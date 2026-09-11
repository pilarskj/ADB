package adb.operators;

import adb.tree.AnnotatedNode;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

// chains parameter scaling and resampling hidden events together
// extends UpDownOperator, modified AnnotatedScale
public class JointEventSampler extends AnnotatedTreeOperator {

    public Input<List<RealParameter>> upInput = new Input<>("up",
            "parameters to scale upwards", new ArrayList<>());

    public Input<List<RealParameter>> downInput = new Input<>("down",
            "parameters to scale downwards", new ArrayList<>());

    public Input<Boolean> elementWiseInput = new Input<>("elementWise", "flag to indicate that the scaling is applied to a random index in multivariate parameters (default false)", false);

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "magnitude factor used for scaling", 0.75);

    public Input<Double> proportionBranchesInput = new Input<>("proportionBranches",
            "proportion of branches to sample events on (default 0.1)", 0.1);

    public Input<Boolean> optimiseInput = new Input<>("optimise",
            "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    double scaleFactor;
    double proportionBranches;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        scaleFactor = scaleFactorInput.get();
        proportionBranches = proportionBranchesInput.get();
        // sanity check
        if (upInput.get().size() + downInput.get().size() == 0) {
            Log.warning.println("WARNING: At least one up or down item must be specified");
        }
    }


    @Override
    public double proposal() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);

        // collect current branch probabilities
        double[] oldBranchProbs = new double[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            oldBranchProbs[i] = getBranchProbability((AnnotatedNode)tree.getNode(i), true);
        }

        // choose scale factor
        double u = Randomizer.nextDouble();
        double f = u * scaleFactor + (1.0-u) / scaleFactor;

        // keep track of Hastings ratio
        double logf = Math.log(f);
        double logHR = -2 * logf;

        // TODO: add option to scale just one element
        // scale parameters upwards
        for (int pidx = 0; pidx < upInput.get().size(); pidx++) {
            RealParameter param = upInput.get().get(pidx);
            for (int i = 0; i < param.getDimension(); i++) {
                double oldValue = param.getValue(i);
                double newValue = oldValue * f;
                if (newValue < param.getLower() || newValue > param.getUpper())
                    return Double.NEGATIVE_INFINITY;

                param.setValue(i, newValue);
                logHR += logf;
            }
        }

        // scale parameters downwards
        for (int pidx = 0; pidx < downInput.get().size(); pidx++) {
            RealParameter param = downInput.get().get(pidx);
            for (int i = 0; i<param.getDimension(); i++) {
                double oldValue = param.getValue(i);
                double newValue = oldValue / f;
                if (newValue < param.getLower() || newValue > param.getUpper())
                    return Double.NEGATIVE_INFINITY;

                param.setValue(i, newValue);
                logHR -= logf;
            }
        }

        // resample events
        double x;
        double newBranchProb;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            x = Randomizer.nextDouble();
            if (x < proportionBranches) {
                newBranchProb = resampleEvents((AnnotatedNode) tree.getNode(i), true);
                logHR += (oldBranchProbs[i] - newBranchProb);
            }
        }
        return logHR;
    }


    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / scaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(double value) {
        scaleFactor = Math.max(Math.min(value, 1.0), 0.0);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }
}
