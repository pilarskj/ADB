import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.parameter.RealParameter;


public class GammaBranchingModel extends SpeciesTreeDistribution {
    final public Input<RealParameter> scaleParameterInput =
            new Input<>("scale", "scale parameter of the Gamma distribution", Input.Validate.REQUIRED);
    final public Input<RealParameter> shapeParameterInput =
            new Input<>("shape", "shape parameter of the Gamma distribution");
    final public Input<RealParameter> rhoParameterInput =
            new Input<>("rho", "sampling probability at the end of the process");
    final public Input<RealParameter> originParameterInput =
            new Input<>("origin", "time of origin of the process");

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
        double rho = rhoParameterInput.get().getValue();
        double t_or = originParameterInput.get().getValue();

        return calculateTreeLogLikelihood(tree, a, b, rho, t_or);
    }

    protected double calculateTreeLogLikelihood(final TreeInterface tree,
                                                final double a, final double b, final double rho,
                                                final double t_or) {

        // stop if tree origin is smaller than root height
        if (tree.getRoot().getHeight() > t_or)
            return Double.NEGATIVE_INFINITY;

        // get branching times from tree
        BranchingTimesList B = BranchingTimes.getBranchingTimes(tree, t_or);
        double[] int_s = B.internalStartTimes;
        double[] int_e = B.internalEndTimes;
        double[] ext_s = B.externalStartTimes;
        double[] ext_e = B.externalEndTimes;
        // assert all(ext_s == 0))

        // set default options for calcLogLikelihood
        int m = (int)Math.pow(2, 14);
        int maxit = 100;

        double logL = GammaLogLikelihood.calcLogLikelihood(rho, a, b, t_or, int_s, int_e, ext_e, m, maxit);
        return logL;
    }

    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation()
                || scaleParameterInput.get().somethingIsDirty()
                || shapeParameterInput.get().somethingIsDirty()
                || rhoParameterInput.get().somethingIsDirty()
                || originParameterInput.get().somethingIsDirty();
    }

    @Override
    public boolean canHandleTipDates() {
        return false;
    }
}
