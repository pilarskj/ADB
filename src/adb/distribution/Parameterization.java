package adb.distribution;

import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;

import java.util.ArrayList;


public class Parameterization extends CalculationNode {

    public Input<Integer> nTypesInput =
            new Input<>("nTypes", "number of types (default 1)", 1);

    public Input<RealParameter> lifetimeParameterInput =
            new Input<>("lifetime", "expected lifetime per type", Input.Validate.REQUIRED);
    public Input<Parameter> shapeParameterInput =
            new Input<>("shape", "shape parameter of the lifetime distribution per type", Input.Validate.REQUIRED);
    public Input<RealParameter> deathParameterInput =
            new Input<>("death", "death probability per type", Input.Validate.REQUIRED);

    public Input<RealParameter> symTransitionsParameterInput =
            new Input<>("symTransitions", "symmetric transition probabilities between types (matrix)", new RealParameter("1.0"));
    public Input<RealParameter> asymTransitionsParameterInput =
            new Input<>("asymTransitions", "asymmetric transition probabilities between types (matrix with zeros on diagonal)", new RealParameter("0.0"));

    public Input<RealParameter> samplingParameterInput =
            new Input<>("sampling", "sampling probability at the end of the process per type", Input.Validate.REQUIRED);

    public Input<Double> originTimeInput =
            new Input<>("originTime", "time since origin of the process"); // TODO: use processLength instead?
    public Input<Integer> originTypeInput =
            new Input<>("originType", "type at origin", 0); // TODO: consider startTypeProbs

    // TODO: allow for custom type values (e.g. "stem", "neuron",...)
    /* public Input<String> typeLabelInput =
            new Input<>("typeLabel", "Attribute key used to specify sample types in tree (default type)", "type");
    public Input<TraitSet> typeSetInput =
            new Input<>("typeSet", "Set specifying sample types."); */


    protected Integer nTypes;
    protected RealParameter lifetime, death, symTransitions, asymTransitions, sampling;
    protected Parameter shape;
    protected Double originTime;
    protected Integer originType;

    boolean dirtyTransitions;
    boolean[][] transitionMatrix, storedTransitionMatrix;
    boolean[][] transitionGraph, storedTransitionGraph;


    @Override
    public void initAndValidate() {

        // get parameters
        nTypes = nTypesInput.get();
        lifetime = lifetimeParameterInput.get();
        shape = shapeParameterInput.get();
        death = deathParameterInput.get();
        symTransitions = symTransitionsParameterInput.get();
        asymTransitions = asymTransitionsParameterInput.get();
        sampling = samplingParameterInput.get();
        originTime = originTimeInput.get();
        originType = originTypeInput.get();

        int nDims = nTypes * nTypes;

        // check length of arrays // TODO: allow single value per input? (same across types)
        if (lifetime.getDimension() != nTypes || shape.getDimension() != nTypes || death.getDimension() != nTypes ||
                symTransitions.getDimension() != nDims || asymTransitions.getDimension() != nDims ||
                sampling.getDimension() != nTypes) {
            throw new IllegalArgumentException("Please check the number of types and parameter dimensions!");
        }

        // check type at origin
        if (originType >= nTypes) {
            throw new IllegalArgumentException("Origin type is not valid.");
        }

        // assert that parameters fall in the proper range
        lifetime.setLower(Math.max(lifetime.getLower(), 0.0));
        shape.setLower(Math.max(((Number)shape.getLower()).doubleValue(), 0));
        death.setBounds(Math.max(death.getLower(), 0.0), Math.min(death.getUpper(), 1.0));
        sampling.setBounds(Math.max(sampling.getLower(), 0.0), Math.min(sampling.getUpper(), 1.0));
        symTransitions.setBounds(Math.max(symTransitions.getLower(), 0.0), Math.min(symTransitions.getUpper(), 1.0));
        asymTransitions.setBounds(Math.max(asymTransitions.getLower(), 0.0), Math.min(asymTransitions.getUpper(), 1.0));

        // check constraints on transition matrices
        for (int i = 0; i < nTypes; i++) {
            if (getAsymTransition(i, i) != 0) {
                throw new IllegalArgumentException("Diagonal entries of the asymmetric transition probability matrix must be 0!");
            }
            double sum = 0;
            for (int j = 0; j < nTypes; j++) {
                sum += getSymTransition(i, j) + getAsymTransition(i, j);
            }
            if (Math.abs(sum - 1.0) > 1e-6) {
                throw new IllegalArgumentException("All transition probabilities per type must sum to 1!");
            }
        }

        // create transition maps
        transitionMatrix = new boolean[nTypes][nTypes];
        storedTransitionMatrix = new boolean[nTypes][nTypes];
        transitionGraph = new boolean[nTypes][nTypes];
        storedTransitionGraph = new boolean[nTypes][nTypes];
        dirtyTransitions = true;
        updateTransitions();
    }

    public int getNTypes() {
        return nTypes;
    }

    public RealParameter getLifetimeParameter() {
        return lifetime;
    }

    public double getLifetime(int i) {
        return lifetime.getArrayValue(i);
    }

    public Parameter getShapeParameter() {
        return shape;
    }

    public boolean shapeIsInteger() {
        return shape instanceof IntegerParameter;
    }

    public double getShape(int i) {
        return shape.getArrayValue(i);
    }

    public double getDeath(int i) {
        return death.getArrayValue(i);
    }

    public double getSampling(int i) {
        return sampling.getArrayValue(i);
    }

    public double getSymTransition(int i, int j) {
        return symTransitions.getArrayValue(i * nTypes + j);
    }

    public double getAsymTransition(int i, int j) {
        return asymTransitions.getArrayValue(i * nTypes + j);
    }

    public Double getOriginTime() {
        return originTime;
    }

    public int getOriginType() {
        return originType;
    }

    public boolean isDirectProgenitor(int i, int j) {
        updateTransitions();
        return transitionMatrix[i][j];
    }

    public boolean isProgenitor(int i, int j) {
        updateTransitions();
        return transitionGraph[i][j];
    }

    public ArrayList<Integer> getDirectProgenitors(int i) {
        updateTransitions();
        ArrayList<Integer> progenitors = new ArrayList<>();
        for (int j = 0; j < nTypes; j++) {
            if (isDirectProgenitor(j, i)) {
                progenitors.add(j);
            }
        }
        return progenitors;
    }

    private void updateTransitions() {
        if (!dirtyTransitions) {
            return;
        }

        // direct transitions (adjacency matrix)
        for (int i = 0; i < nTypes; i++) {
            for (int j = 0; j < nTypes; j++) {
                if (getSymTransition(i, j) > 0 || getAsymTransition(i, j) > 0) {
                    transitionMatrix[i][j] = true;
                }
            }
        }
        // indirect transitions (transitive closure) -- Warshall's Algorithm, cf. https://cs.winona.edu/lin/cs440/ch08-2.pdf
        // deep copy for intialisation
        for (int i = 0; i < nTypes; i++) {
            System.arraycopy(transitionMatrix[i], 0, transitionGraph[i], 0, nTypes);
        }
        // iteration
        for (int k = 0; k < nTypes; k++) {
            for (int i = 0; i < nTypes; i++) {
                for (int j = 0; j < nTypes; j++) {
                    transitionGraph[i][j] = transitionGraph[i][j] || (transitionGraph[i][k] && transitionGraph[k][j]);
                }
            }
        }

        dirtyTransitions = false;
    }

    @Override
    public boolean requiresRecalculation() {
        if (symTransitions.somethingIsDirty() || asymTransitions.somethingIsDirty()) {
            // TODO: monitor only binary changes (0 or >0), not probability updates
            dirtyTransitions = true;
        }

        return true;
    }


    @Override
    protected void store() {
        for (int i = 0; i < nTypes; i++) {
            System.arraycopy(transitionMatrix[i], 0, storedTransitionMatrix[i], 0, nTypes);
            System.arraycopy(transitionGraph[i], 0, storedTransitionGraph[i], 0, nTypes);
        }
        super.store();
    }


    @Override
    protected void restore() {
        boolean[][] tmp;

        tmp = transitionMatrix;
        transitionMatrix = storedTransitionMatrix;
        storedTransitionMatrix = tmp;

        tmp = transitionGraph;
        transitionGraph = storedTransitionGraph;
        storedTransitionGraph = tmp;

        super.restore();
    }

}
