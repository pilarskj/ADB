package adb.operators;

import adb.distribution.Parameterization;
import adb.tree.AnnotatedNode;
import beast.base.core.Input;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

// chains lifetime scaling and resampling hidden events together
// currently partly redundant with AnnotatedScale
public class JointEventSampler extends AnnotatedTreeOperator {

    public Input<RealParameter> realParameterInput =
            new Input<>("realParameter", "parameter to be moved (either shape or lifetime)");

    public Input<IntegerParameter> intParameterInput =
            new Input<>("intParameter", "parameter to be moved (either shape or lifetime)", Input.Validate.XOR, realParameterInput);

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "for scaling real parameter, is restricted to the range [1/scaleFactor, scaleFactor]", 0.75);

    public Input<Integer> windowSizeInput = new Input<>("windowSize",
            "the size of the window both up and down", 1);

    public Input<Boolean> drawEventCountInput = new Input<>("drawEventCount",
            "draw or fix the number of events along internal branches (default true)", true);

    public Input<Boolean> optimiseInput = new Input<>("optimise",
            "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    double scaleFactor;
    int windowSize;
    boolean drawEventCount;

    boolean isReal;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        scaleFactor = scaleFactorInput.get();
        windowSize = windowSizeInput.get();
        drawEventCount = drawEventCountInput.get();
        isReal = realParameterInput.get() != null;
        // TODO: some checks for compatibility?
    }


    @Override
    public double proposal() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);

        double logHR = 0.0;
        double newBranchProb;

        // TODO: first check validity of parameter move?
        // TODO: keep track of indicators for multi-type case

        // collect current branch probabilities
        double[] oldBranchProbs = new double[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            oldBranchProbs[i] = getBranchProbability((AnnotatedNode)tree.getNode(i), drawEventCount);
        }

        if (isReal) {
            RealParameter realParameter = (RealParameter) InputUtil.get(realParameterInput,this);

            // choose scale factor
            double u = Randomizer.nextDouble();
            double f = u * scaleFactor + (1.0-u) / scaleFactor;

            // keep track of Hastings ratio
            double logf = Math.log(f);
            logHR += -2 * logf;

            // scale
            double oldValue = realParameter.getValue(0);
            double newValue = oldValue * f;

            if (newValue < realParameter.getLower() || newValue > realParameter.getUpper())
                return Double.NEGATIVE_INFINITY;

            realParameter.setValue(0, newValue);
            logHR += logf;

        } else {
            IntegerParameter intParameter = (IntegerParameter) InputUtil.get(intParameterInput,this);

            // choose magnitude of move
            // int i = Randomizer.nextInt(param.getDimension());
            int oldValue = intParameter.getValue(0);
            int newValue = oldValue + Randomizer.nextInt(2 * windowSize + 1) - windowSize;

            if (newValue < intParameter.getLower() || newValue > intParameter.getUpper()) {
                return Double.NEGATIVE_INFINITY;
            }
            if (newValue == oldValue) {
                return Double.NEGATIVE_INFINITY;
            }

            intParameter.setValue(0, newValue);
        }

        // resample events
        for (int i = 0; i < tree.getNodeCount(); i++) {
            newBranchProb = resampleEvents((AnnotatedNode)tree.getNode(i), drawEventCount);
            logHR += (oldBranchProbs[i] - newBranchProb);
        }
        // for testing
        // Tree flatTree = ((AnnotatedTree)tree).convertAnnotatedTree(false);
        // System.out.println(flatTree.getRoot().toNewick());
        return logHR;
    }

    // TODO: add optimisation of scaling factor?
}
