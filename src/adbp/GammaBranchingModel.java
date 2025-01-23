package adbp;

import beast.base.core.*;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;

import static org.apache.commons.math3.special.Gamma.logGamma;

@Description("This model implements an Age-Dependent Branching Process " +
        "with lifetimes distributed according to a Gamma distribution with integer shape parameter (Erlang distribution), " +
        "a death probability and extant sampling.")
public class GammaBranchingModel extends SpeciesTreeDistribution {

    // parametrization
    final public Input<RealParameter> lifetimeParameterInput =
            new Input<>("lifetime", "expected lifetime", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> shapeIntegerParameterInput =
            new Input<>("shapeInteger", "integer shape parameter of the Gamma distribution");
    final public Input<RealParameter> shapeRealParameterInput =
            new Input<>("shapeReal", "real shape parameter of the Gamma distribution");
    final public Input<RealParameter> deathParameterInput =
            new Input<>("deathprob", "probability of death at branching times (default 0)", new RealParameter("0.0"));
    final public Input<RealParameter> rhoParameterInput =
            new Input<>("rho", "sampling probability at the end of the process (default 1)", new RealParameter("1.0"));
    final public Input<RealParameter> originParameterInput =
            new Input<>("origin", "time of origin of the process", Input.Validate.REQUIRED);

    // options
    public Input<Integer> maxIterationsInput =
            new Input<>("maxIterations", "maximum number of iterations for numerical integration",100);
    public Input<Double> tolerancePInput =
            new Input<>("toleranceP", "tolerance for numerical integration of P0 and P1",1e-12);
    public Input<Double> toleranceBInput =
            new Input<>("toleranceB", "tolerance for numerical integration of B (branch probabilities)",1e-6);
    public Input<Integer> stepSizePInput =
            new Input<>("stepSizeP", "number of time steps for FFT for P0 and P1", (int)Math.pow(2, 14));
    public Input<Integer> stepSizeBInput =
            new Input<>("stepSizeB", "number of time steps for FFT for B", (int)Math.pow(2, 12));
    public Input<Boolean> approxInput =
            new Input<>("approx", "approximate branch probabilities (default true)", true);


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // make sure that all tips are at the same height
        TreeInterface tree = treeInput.get();
        if (!TreeUtils.isUltrametric(tree)) {
            Log.warning.println("WARNING: This model (tree prior) cannot handle dated tips.");
        }

        // check type of shape
        if (shapeIntegerParameterInput.get() == null && shapeRealParameterInput.get() == null) {
            throw new IllegalArgumentException("Please specify a shape parameter");
        } else if (shapeIntegerParameterInput.get() != null && shapeRealParameterInput.get() != null) {
            throw new IllegalArgumentException("Please specify either an integer or real shape parameter, not both");
        } else if (shapeRealParameterInput.get() != null) {
            if (approxInput.get()) {
                throw new IllegalArgumentException("The approximation only works for an integer shape parameter!");
            }
        }

        // check that step sizes are power of 2 (required for FFT)
        if (!isPowerOfTwo(stepSizePInput.get())) {
            throw new IllegalArgumentException("stepSizeP must be a power of 2");
        }

        if (!isPowerOfTwo(stepSizeBInput.get())) {
            throw new IllegalArgumentException("stepSizeB must be a power of 2");
        }
    }


    // Helper method to check if number is a power if two
    private boolean isPowerOfTwo(int x) {
        return (x > 0) && ((x & (x - 1)) == 0); // https://stackoverflow.com/questions/600293/how-to-check-if-a-number-is-a-power-of-2
    }


    @Override
    public double calculateTreeLogLikelihood(final TreeInterface tree) {

        // parameters
        double C = lifetimeParameterInput.get().getValue();
        Number b = null;
        if (shapeIntegerParameterInput.get() != null) {
            b = shapeIntegerParameterInput.get().getValue();
        } else if (shapeRealParameterInput.get() != null) {
            b = shapeRealParameterInput.get().getValue();
        }
        double d = deathParameterInput.get().getValue();
        double rho = rhoParameterInput.get().getValue();
        double origin = originParameterInput.get().getValue();

        // options
        int maxIt = maxIterationsInput.get();
        double tolP = tolerancePInput.get();
        double tolB = toleranceBInput.get();
        int mP = stepSizePInput.get();
        int mB = stepSizeBInput.get();
        boolean approx = approxInput.get();

        return calculateTreeLogLikelihood(tree, C, b, d, rho, origin, maxIt, tolP, tolB, mP, mB, approx);
    }


    protected double calculateTreeLogLikelihood(final TreeInterface tree,
                                                final double C, final Number b, final double d, final double rho, final double origin,
                                                final int maxIt, final double tolP, final double tolB, final int mP, final int mB, final boolean approx) {

        // stop if tree origin is smaller than root height
        if (tree.getRoot().getHeight() > origin)
            return Double.NEGATIVE_INFINITY;

        // get branching times from tree
        // create empty arrays to store start and end times of branches
        double[] intS = new double[0];
        double[] intE = new double[0];
        double[] extE = new double[0];

        // traverse all nodes and compute start/end times of all branches (backwards in time)
        // alternatively, use Branch/ BranchList classes
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);

            // start time is the node height (or divergence time of the node)
            double startTime = node.getHeight();

            double endTime;
            if (node.isRoot()) {
                // end time is the tree origin
                endTime = origin;
            } else {
                // end time is the parent node height (i.e., the divergence time of the parent)
                endTime = node.getParent().getHeight();
            }

            // if node is a tip, add times to external branches, otherwise, add to internal branches
            if (node.isLeaf()) {
                assert startTime == 0; // tips should have height 0
                extE = Arrays.copyOf(extE, extE.length + 1);
                extE[extE.length - 1] = endTime;
            } else {
                intS = Arrays.copyOf(intS, intS.length + 1);
                intE = Arrays.copyOf(intE, intE.length + 1);
                intS[intS.length - 1] = startTime;
                intE[intE.length - 1] = endTime;
            }
        }

        // check node numbers
        int nTips = tree.getLeafNodeCount();
        assert extE.length == nTips;
        assert intS.length == tree.getInternalNodeCount();
        assert intE.length == tree.getInternalNodeCount();

        // calculate tree factor
        double treeFactor = (nTips - 1) * Math.log(2) - logGamma(nTips + 1); // 2^(n-1)/n!

        // get scale parameter of the Gamma distribution
        double a = C / b.doubleValue();

        // calculate LogLikelihood
        double logL = GammaLogLikelihood.calcLogLikelihood(a, b, d, rho, origin,
                intS, intE, extE, maxIt, tolP, tolB, mP, mB, approx);

        return treeFactor + logL;
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
