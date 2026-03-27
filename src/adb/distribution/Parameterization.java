package adb.distribution;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;

public class Parameterization extends CalculationNode {

    public Input<Integer> nTypesInput =
            new Input<>("nTypes", "number of types (default 1)", 1);

    public Input<RealParameter> lifetimeParameterInput =
            new Input<>("lifetime", "expected lifetime per type", Input.Validate.REQUIRED);
    public Input<Parameter> shapeParameterInput =
            new Input<>("shape", "shape parameter of the lifetime distribution per type", Input.Validate.REQUIRED);
    public Input<RealParameter> deathParameterInput =
            new Input<>("death", "death probability per type", Input.Validate.REQUIRED);

    public Input<RealParameter> transitionsParameterInput =
            new Input<>("transitions", "transition probabilities between types (matrix)");
    public Input<RealParameter> sTransitionsParameterInput =
            new Input<>("sTransitions", "symmetric transition probabilities between types (matrix)");
    public Input<RealParameter> asTransitionsParameterInput =
            new Input<>("asTransitions", "asymmetric transition probabilities between types (matrix with zeros on diagonal)");

    public Input<RealParameter> samplingParameterInput =
            new Input<>("sampling", "sampling probability at the end of the process per type", Input.Validate.REQUIRED);

    public Input<Double> originTimeInput =
            new Input<>("originTime", "time since origin of the process"); // TODO: use processLength instead?
    public Input<Integer> originTypeInput =
            new Input<>("originType", "type at origin"); // TODO: consider startTypeProbs

    // TODO: allow for custom type values (e.g. "stem", "neuron",...)
    /* public Input<String> typeLabelInput =
            new Input<>("typeLabel", "Attribute key used to specify sample types in tree (default type)", "type");
    public Input<TraitSet> typeSetInput =
            new Input<>("typeSet", "Set specifying sample types."); */


    protected int nTypes;
    protected double[] lifetime, shape, death, sampling; // TODO: do I need Number[] or int[] shape in approximate likelihood calculation?
    protected double[][] sTransitions, asTransitions;
    protected Double originTime;
    protected int originType;

    protected double[] storedLifetime, storedShape, storedDeath, storedSampling;
    protected double[][] storedSTransitions, storedAsTransitions;


    @Override
    public void initAndValidate() {

        // get parameters
        nTypes = nTypesInput.get();
        lifetime = lifetimeParameterInput.get().getDoubleValues();
        shape = shapeParameterInput.get().getDoubleValues();
        death = deathParameterInput.get().getDoubleValues();
        sampling = samplingParameterInput.get().getDoubleValues();
        originTime = originTimeInput.get();

        // check length of arrays // TODO: allow single value per input? (same across types)
        if (lifetime.length != nTypes || shape.length != nTypes || death.length != nTypes || sampling.length != nTypes) {
            throw new IllegalArgumentException("Please check the number of types and parameter dimensions!");
        }

        // TODO: are these loops performed when any of the parameter changes or only at the start of MCMC?
        if (nTypes == 1) {
            if (transitionsParameterInput.get() != null ||
                    sTransitionsParameterInput.get() != null || asTransitionsParameterInput.get() != null ||
                    originTypeInput.get() != null) {
                Log.warning("WARNING: Initializing single-type model, provided transition probabilities and type at origin will be ignored.");
            }
            sTransitions = new double[][]{{1.0}};
            asTransitions = new double[][]{{0.0}};
            originType = 0;

        } else {
            originType = originTypeInput.get();
            if (originType >= nTypes) {
                throw new IllegalArgumentException("Origin type is not valid.");
            }

            int nDims = nTypes * nTypes;
            sTransitions = new double[nTypes][nTypes];
            asTransitions = new double[nTypes][nTypes];

            if (transitionsParameterInput.get() == null && sTransitionsParameterInput.get() != null && asTransitionsParameterInput.get() != null) {
                double[] sTransitionsArray = sTransitionsParameterInput.get().getDoubleValues();
                double[] asTransitionsArray = asTransitionsParameterInput.get().getDoubleValues();
                if (sTransitionsArray.length != nDims || asTransitionsArray.length != nDims) {
                    throw new IllegalArgumentException("Please check the dimension of transition matrices!");
                }
                // fill matrices and assert that constraints are met
                for (int i = 0; i < nTypes; i++) {
                    double sum = 0;
                    for (int j = 0; j < nTypes; j++) {
                        sTransitions[i][j] = sTransitionsArray[i * nTypes + j];
                        asTransitions[i][j] = asTransitionsArray[i * nTypes + j];
                        sum += sTransitions[i][j] + asTransitions[i][j];
                    }
                    if (asTransitions[i][i] != 0) {
                        throw new IllegalArgumentException("Diagonal entries of the asymmetric transition probability matrix must be 0!");
                    }
                    if (Math.abs(sum - 1.0) > 1e-6) {
                        throw new IllegalArgumentException("All transition probabilities per type must sum to 1!");
                    }
                }

            } else if (transitionsParameterInput.get() != null && sTransitionsParameterInput.get() == null && asTransitionsParameterInput.get() == null) {
                double[] transitionsArray = transitionsParameterInput.get().getDoubleValues();
                if (transitionsArray.length != nDims) {
                    throw new IllegalArgumentException("Please check the transition matrix dimension!");
                }
                // fill matrices and assert that constraints are met
                for (int i = 0; i < nTypes; i++) {
                    double sum = 0;
                    for (int j = 0; j < nTypes; j++) {
                        double value = transitionsArray[i * nTypes + j];
                        if (j == i) {
                            sTransitions[i][j] = value;
                        } else {
                            asTransitions[i][j] = value;
                        }
                        sum += value;
                    }
                    if (Math.abs(sum - 1.0) > 1e-6) {
                        throw new IllegalArgumentException("All transition probabilities per type must sum to 1!");
                    }
                }

            } else {
                throw new IllegalArgumentException("Either provide one simplified transition matrix, " +
                        "or both symmetric and asymmetric transition matrices.");
            }
        }

        // TODO: how are matrices updated when parameters change?
    }
}
