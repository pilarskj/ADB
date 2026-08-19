package adb.operators;

import adb.tree.AnnotatedNode;
import beast.base.core.Input;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

// chains lifetime scaling and resampling hidden events together
// currently partly redundant with AnnotatedScale
public class JointEventSampler extends AnnotatedTreeOperator {

    public Input<List<RealParameter>> parametersInput =
            new Input<>("parameter",
                    "scale this scalar parameter",
                    new ArrayList<RealParameter>());

    public Input<List<RealParameter>> parametersInverseInput =
            new Input<>("parameterInverse",
                    "scale this scalar parameter inversely",
                    new ArrayList<RealParameter>());

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "for scaling real parameter, is restricted to the range [1/scaleFactor, scaleFactor]", 0.75);

    public Input<Boolean> optimiseInput = new Input<>("optimise",
            "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    double scaleFactor;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        scaleFactor = scaleFactorInput.get();
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

        // scale parameters
        for (int pidx=0; pidx<parametersInput.get().size(); pidx++) {
            RealParameter param = parametersInput.get().get(pidx);
            for (int i=0; i<param.getDimension(); i++) {
                double oldValue = param.getValue(i);
                double newValue = oldValue * f;
                if (newValue < param.getLower() || newValue > param.getUpper())
                    return Double.NEGATIVE_INFINITY;

                param.setValue(i, newValue);
                logHR += logf;
            }
        }

        // scale parameters inversely
        for (int pidx=0; pidx<parametersInverseInput.get().size(); pidx++) {
            RealParameter param = parametersInverseInput.get().get(pidx);
            for (int i=0; i<param.getDimension(); i++) {
                double oldValue = param.getValue(i);
                double newValue = oldValue/f;
                if (newValue < param.getLower() || newValue > param.getUpper())
                    return Double.NEGATIVE_INFINITY;

                param.setValue(i, newValue);
                logHR -= logf;
            }
        }

        // resample events
        double newBranchProb;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            newBranchProb = resampleEvents((AnnotatedNode)tree.getNode(i), true);
            logHR += (oldBranchProbs[i] - newBranchProb);
        }
        return logHR;
    }

    // TODO: add optimisation of scaling factor?
}
