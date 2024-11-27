package adbp;

import beast.base.core.*;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;

@Description("This model implements an Age-Dependent Branching Process " +
        "with lifetimes distributed according to a Gamma distribution with integer shape parameter (Erlang distribution), " +
        "a death probability and extant sampling.")
public class GammaBranchingModel extends SpeciesTreeDistribution {

    // parametrization
    final public Input<RealParameter> scaleParameterInput =
            new Input<>("scale", "scale parameter of the Gamma distribution", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> shapeParameterInput =
            new Input<>("shape", "shape parameter of the Gamma distribution", Input.Validate.REQUIRED);
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
        // Ugne: mP and mB should always be power of 2 (per your comment elsewhere).
        // I would add a check for this here. For example:
        // if (!isPowerOfTwo(stepSizePInput.get())) {
        //        throw new IllegalArgumentException("stepSizeP must be a power of 2");
        //    }
        //
        // if (!isPowerOfTwo(stepSizeBInput.get())) {
        //        throw new IllegalArgumentException("stepSizeB must be a power of 2");
        //    }
        // where
        // private boolean isPowerOfTwo(int x) {
        //      return (x > 0) && ((x & (x - 1)) == 0); // https://stackoverflow.com/questions/600293/how-to-check-if-a-number-is-a-power-of-2
        // }
    }


    @Override
    public double calculateTreeLogLikelihood(final TreeInterface tree) {

        // parameters
        double a = scaleParameterInput.get().getValue();
        int b = shapeParameterInput.get().getValue();
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

        return calculateTreeLogLikelihood(tree, a, b, d, rho, origin, maxIt, tolP, tolB, mP, mB, approx);
    }


    protected double calculateTreeLogLikelihood(final TreeInterface tree,
                                                final double a, final int b, final double d, final double rho, final double origin,
                                                final int maxIt, final double tolP, final double tolB, final int mP, final int mB, final boolean approx) {

        // stop if tree origin is smaller than root height
        if (tree.getRoot().getHeight() > origin)
            return Double.NEGATIVE_INFINITY;

        // get branching times from tree
        // create empty arrays to store start and end times of branches
        double[] int_s = new double[0];
        double[] int_e = new double[0];
        double[] ext_e = new double[0];

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
                ext_e = Arrays.copyOf(ext_e, ext_e.length + 1);
                ext_e[ext_e.length - 1] = endTime;
            } else {
                int_s = Arrays.copyOf(int_s, int_s.length + 1);
                int_e = Arrays.copyOf(int_e, int_e.length + 1);
                int_s[int_s.length - 1] = startTime;
                int_e[int_e.length - 1] = endTime;
            }
        }

        // calculate LogLikelihood
        double logL = GammaLogLikelihood.calcLogLikelihood(a, b, d, rho, origin,
                int_s, int_e, ext_e, maxIt, tolP, tolB, mP, mB, approx);

        return logL;
    }


    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation()
                || scaleParameterInput.get().somethingIsDirty()
                || shapeParameterInput.get().somethingIsDirty()
                || deathParameterInput.get().somethingIsDirty()
                || rhoParameterInput.get().somethingIsDirty()
                || originParameterInput.get().somethingIsDirty();
    }


    @Override
    public boolean canHandleTipDates() {
        return false;
    }
}
