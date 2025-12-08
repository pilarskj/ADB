package adb;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;

import static adb.MTLogLikelihood.calcMTLogLikelihood;
import static org.apache.commons.math3.special.Gamma.logGamma;


@Description("This class implements a multi-type Age-Dependent Branching (mtADB) model " +
        "with type-specific Gamma-distributed lifetimes and death probabilities," +
        "symmetric and asymmetric type transitions, and extant sampling.")
public class MTBranchingModel extends SpeciesTreeDistribution { // TODO: change name to mtADBModel? + Remove single-type classes (ultimate goal: one model for all settings)

    // parametrization
    public Input<RealParameter> lifetimeParameterInput =
            new Input<>("lifetime", "mean lifetime per type (array)", Input.Validate.REQUIRED);
    public Input<RealParameter> shapeParameterInput =
            new Input<>("shape", "shape parameters of the lifetime distribution per type (array)", Input.Validate.REQUIRED);
    public Input<RealParameter> deathParameterInput =
            new Input<>("death", "death probability per type (array)", Input.Validate.REQUIRED);
    public Input<RealParameter> sTransitionsParameterInput =
            new Input<>("sTransitions", "symmetric transition probabilities between types (matrix)", Input.Validate.REQUIRED);
    public Input<RealParameter> asTransitionsParameterInput =
            new Input<>("asTransitions", "asymmetric transition probabilities between types (matrix with zeros on diagonal)", Input.Validate.REQUIRED);
    public Input<RealParameter> rhoParameterInput =
            new Input<>("rho", "sampling probability at the end of the process (default 1)", new RealParameter("1.0")); // TODO: consider type-specific sampling

    public Input<Double> originTimeInput =
            new Input<>("originTime", "time of origin of the process"); // TODO: should origin be a variable or fixed parameter?
    public Input<Integer> originTypeInput =
            new Input<>("originType", "type at origin", 0); // TODO: consider startTypeProbs instead
    // TODO: allow for custom type values (e.g. "stem", "neuron",...)
    /* public Input<String> typeLabelInput =
            new Input<>("typeLabel", "Attribute key used to specify sample types in tree (default type)", "type");
    public Input<TraitSet> typeSetInput =
            new Input<>("typeSet", "Set specifying sample types."); */

    // computing options
    public Input<Integer> maxIterationsInput =
            new Input<>("maxIterations", "maximum number of iterations for numerical integration",100);
    public Input<Double> tolerancePInput =
            new Input<>("toleranceP", "tolerance for numerical integration of P0 and P1",1e-12);
    public Input<Double> toleranceBInput =
            new Input<>("toleranceB", "tolerance for numerical integration of B (branch probabilities)",1e-6);
    public Input<Integer> stepSizePInput =
            new Input<>("stepSizeP", "number of time steps for FFT for P0 and P1", (int)Math.pow(2, 12));
    public Input<Integer> stepSizeBInput =
            new Input<>("stepSizeB", "number of time steps for FFT for B", (int)Math.pow(2, 10));
    public Input<Boolean> conditionOnRootInput =
            new Input<>("conditionOnRoot", "condition on the root height otherwise on origin (default false)", false); // TODO: implement! (consider unknown type at root)

    // define objects
    private TreeInterface tree;
    private BranchList branches;
    private int ntips;
    private int ntypes;
    private double[] lifetime;
    private double[] shape;
    private double[] death;
    private double[][] sTransitions;
    private double[][] asTransitions;
    private double rho;

    private Double originTime;
    private int originType;

    private int maxIt;
    private double tolP;
    private double tolB;
    private int mP;
    private int mB;
    private boolean conditionOnRoot;


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        tree = treeInput.get();
        ntips = tree.getLeafNodeCount();

        // make sure that all tips are at the same height
        if (!TreeUtils.isUltrametric(tree)) {
            Log.warning.println("WARNING: This model (tree prior) cannot handle dated tips.");
        }

        originTime = originTimeInput.get();
        conditionOnRoot = conditionOnRootInput.get();

        // check that origin is given if not conditioning on the root
        if (originTime == null && !conditionOnRoot) {
            throw new IllegalArgumentException("Please provide the time of origin of the process, or use conditionOnRoot!");
        }

        // get parameters
        lifetime = lifetimeParameterInput.get().getDoubleValues();
        shape = shapeParameterInput.get().getDoubleValues();
        death = deathParameterInput.get().getDoubleValues();
        double[] sTransitionsArray = sTransitionsParameterInput.get().getDoubleValues();
        double[] asTransitionsArray = asTransitionsParameterInput.get().getDoubleValues();
        rho = rhoParameterInput.get().getValue();

        ntypes = lifetime.length;
        int ndims = ntypes * ntypes;

        // check length of arrays
        if (shape.length != ntypes || death.length != ntypes || sTransitionsArray.length != ndims || asTransitionsArray.length != ndims) {
            throw new IllegalArgumentException("Please check the dimensions of parameters!");
        }

        // fill matrices and assert that constraints are met
        sTransitions = new double[ntypes][ntypes];
        asTransitions = new double[ntypes][ntypes];
        for (int i = 0; i < ntypes; i++) {
            double sum = 0;
            for (int j = 0; j < ntypes; j++) {
                sTransitions[i][j] = sTransitionsArray[i * ntypes + j];
                asTransitions[i][j] = asTransitionsArray[i * ntypes + j];
                sum += sTransitions[i][j] + asTransitions[i][j];
            }
            if (asTransitions[i][i] != 0) {
                throw new IllegalArgumentException("Diagonal entries of the asymmetric transition probability matrix must be 0!");
            }
            if (Math.abs(sum - 1.0) > 1e-6) {
                throw new IllegalArgumentException("All transition probabilities per type must sum to 1!");
            }
        }

        // check if type at origin valid // TODO: currently assuming that type values = type indices, i.e. 0,1,...,n-1 - generalize!
        originType = originTypeInput.get();
        if (originType >= ntypes) {
            throw new IllegalArgumentException("Type at origin must be a valid type!");
        }

        // get options
        maxIt = maxIterationsInput.get();
        tolP = tolerancePInput.get();
        tolB = toleranceBInput.get();
        mP = stepSizePInput.get();
        mB = stepSizeBInput.get();

        // check that step sizes are power of 2 (required for FFT)
        if (!isPowerOfTwo(mP) || !isPowerOfTwo(mB)) {
            throw new IllegalArgumentException("Step sizes for FFT must be a power of 2!");
        }
    }


    // Helper method to check if number is a power if two
    private boolean isPowerOfTwo(int x) {
        return (x > 0) && ((x & (x - 1)) == 0); // https://stackoverflow.com/questions/600293/how-to-check-if-a-number-is-a-power-of-2
    }


    @Override
    public double calculateTreeLogLikelihood(final TreeInterface tree) {

        Node root = tree.getRoot();
        double rootHeight = root.getHeight();

        // stop if tree origin is smaller than root height
        if (originTime != null) {
            if (rootHeight >= originTime) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // calculate tree factor
        // double treeFactor = (ntips - 1) * Math.log(2) - logGamma(ntips + 1); // 2^(n-1)/n! (currently in LogL function)

        // get scale parameters (theta)
        double[] scale = new double[ntypes];
        for (int i = 0; i < ntypes; i++) {
            scale[i] = lifetime[i] / shape [i];
        }

        // get list of branches in the tree (annotated with types)
        branches = new BranchList(tree, originTime);

        // calculate likelihood
        double logL = MTLogLikelihood.calcMTLogLikelihood(scale, shape, death, rho, sTransitions, asTransitions, originTime, originType, branches,
                maxIt, tolP, tolB, mP, mB);

        return logL;
    }


    @Override
    protected boolean requiresRecalculation() {
        return true;
    }


    @Override
    public boolean canHandleTipDates() {
        return false;
    }
}
