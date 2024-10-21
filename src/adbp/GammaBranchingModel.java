package adbp;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;


public class GammaBranchingModel extends SpeciesTreeDistribution {
    final public Input<RealParameter> scaleParameterInput =
            new Input<>("scale", "scale parameter of the Gamma distribution", Input.Validate.REQUIRED);
    final public Input<RealParameter> shapeParameterInput =
            new Input<>("shape", "shape parameter of the Gamma distribution");
    final public Input<RealParameter> deathParameterInput =
            new Input<>("deathprob", "probability of death at branching times");
    final public Input<RealParameter> rhoParameterInput =
            new Input<>("rho", "sampling probability at the end of the process");
    final public Input<RealParameter> originParameterInput =
            new Input<>("origin", "time of origin of the process");
    /*
    // for testing time steps and efficiency
    final public Input<IntegerParameter> mPInput =
            new Input<>("mP", "");
    final public Input<IntegerParameter> mBInput =
            new Input<>("mB", "");
    */

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // make sure that all tips are at the same height
        TreeInterface tree = treeInput.get();
        if (!TreeUtils.isUltrametric(tree)) {
            Log.warning.println("WARNING: This model (tree prior) cannot handle dated tips.");
        }
    }

    @Override
    public double calculateTreeLogLikelihood(final TreeInterface tree) {

        double a = scaleParameterInput.get().getValue();
        double b = shapeParameterInput.get().getValue();
        double d = deathParameterInput.get().getValue();
        double rho = rhoParameterInput.get().getValue();
        double t_or = originParameterInput.get().getValue();
        //int mP = mPInput.get().getValue();
        //int mB = mBInput.get().getValue();

        return calculateTreeLogLikelihood(tree, a, b, d, rho, t_or); //mP, mB
    }

    protected double calculateTreeLogLikelihood(final TreeInterface tree,
                                                final double a, final double b, final double d, final double rho,
                                                final double t_or) {
                                                //final int mP, final int mB

        // stop if tree origin is smaller than root height
        if (tree.getRoot().getHeight() > t_or)
            return Double.NEGATIVE_INFINITY;

        // set default options for calcLogLikelihood
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);
        int maxit = 100;

        // get branching times from tree
        // create empty arrays to store start and end times of branches
        double[] int_s = new double[0];
        double[] int_e = new double[0];
        double[] ext_e = new double[0];

        // traverse all nodes and compute start/end times of all branches (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);

            // start time is the node height (or divergence time of the node)
            double startTime = node.getHeight();

            double endTime;
            if (node.isRoot()) {
                // end time is the tree origin
                endTime = t_or;
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
        double logL = GammaLogLikelihood.calcLogLikelihood(rho, a, b, d, t_or, int_s, int_e, ext_e, mP, mB, maxit);
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
