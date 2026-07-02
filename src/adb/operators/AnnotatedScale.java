package adb.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedNode.EventNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

// adapted from https://github.com/tgvaughan/MultiTypeTree/blob/master/src/multitypetree/operators/MultiTypeTreeScale.java
// assumes ultrametric trees (leaf height = 0.0!)
@Description("Scale operator for annotated trees. Also allows additional "
        + "scalar parameters to be rescaled (either forward or inversely) "
        + "at the same time.")
public class AnnotatedScale extends AnnotatedTreeOperator {
    
    public Input<List<RealParameter>> parametersInput =
            new Input<>("parameter",
            "scale this scalar parameter by the same amount as tree",
            new ArrayList<RealParameter>());
    
    public Input<List<BooleanParameter>> indicatorsInput =
            new Input<>("indicator",
            "if provided, used to specify a subset of parameter elements to scale",
            new ArrayList<BooleanParameter>());
    
    public Input<List<RealParameter>> parametersInverseInput =
            new Input<>("parameterInverse",
            "scale this scalar parameter inversely",
            new ArrayList<RealParameter>());
    
    public Input<List<BooleanParameter>> indicatorsInverseInput =
            new Input<>("indicatorInverse",
            "if provided, used to specify a subset of parameter elements to scale inversely",
            new ArrayList<BooleanParameter>());

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "scaling is restricted to the range [1/scaleFactor, scaleFactor]", 0.75);

    final public Input<Boolean> optimiseInput = new Input<>("optimise",
            "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    public Input<Boolean> drawEventsInput = new Input<>("drawEvents",
            "re-draw events along branches (default true)", true);

    boolean drawEvents;
    double scaleFactor;
    boolean indicatorsUsed, indicatorsInverseUsed;
    
    @Override
    public void initAndValidate() {

        super.initAndValidate();
        scaleFactor = scaleFactorInput.get();
        drawEvents = drawEventsInput.get();

        // TODO: keep or remove?
        if (indicatorsInput.get().size()>0) {
            if (indicatorsInput.get().size() != parametersInput.get().size())
                throw new IllegalArgumentException("If an indicator element "
                        + "exists, the number of such elements must equal "
                        + "the number of parameter elements.");
            
            for (int pidx=0; pidx<parametersInput.get().size(); pidx++) {
                if (parametersInput.get().get(pidx).getDimension() != 
                        indicatorsInput.get().get(pidx).getDimension()) {
                    throw new IllegalArgumentException("The number of boolean "
                            + "values in indicator element "
                            + String.valueOf(pidx+1)
                            + " doesn't match the dimension of the "
                            + "corresponding parameter element.");
                }
            }
            indicatorsUsed = true;
        } else
            indicatorsUsed = false;
        
        if (indicatorsInverseInput.get().size()>0) {
            if (indicatorsInverseInput.get().size() != parametersInverseInput.get().size())
                throw new IllegalArgumentException("If an indicatorInverse element "
                        + "exists, the number of such elements must equal "
                        + "the number of parameterInverse elements.");
            
            for (int pidx=0; pidx<parametersInverseInput.get().size(); pidx++) {
                if (parametersInverseInput.get().get(pidx).getDimension() != 
                        indicatorsInverseInput.get().get(pidx).getDimension()) {
                    throw new IllegalArgumentException("The number of boolean "
                            + "values in indicatorInverse element "
                            + String.valueOf(pidx+1)
                            + " doesn't match the dimension of the "
                            + "corresponding parameterInverse element.");
                }
            }
            indicatorsInverseUsed = true;
        } else
            indicatorsInverseUsed = false;
    }


    @Override
    public double proposal() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);

        // choose scale factor
        double u = Randomizer.nextDouble();
        double f = u * scaleFactor + (1.0-u) / scaleFactor;

        // keep track of Hastings ratio
        double logf = Math.log(f);
        double logHR = -2 * logf;

        // scaling cannot reach above origin
        if (parameterization.getOriginTime() != null) {
            AnnotatedNode root = (AnnotatedNode)tree.getRoot();
            if (drawEvents) {
                if (tree.getRoot().getHeight() * f > parameterization.getOriginTime()) {
                    return Double.NEGATIVE_INFINITY;
                }
            } else {
                if (root.getEvent(root.getEventCount()-1).getHeight() * f > parameterization.getOriginTime()) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }

        // collect current branch probabilities (in case of re-drawing events)
        double[] oldBranchProbs = new double[tree.getNodeCount()];
        if (drawEvents) {
            for (int i = 0; i < tree.getNodeCount(); i++) {
                oldBranchProbs[i] = getBranchProbability((AnnotatedNode)tree.getNode(i), true);
            }
        }

        // scale parameters
        for (int pidx=0; pidx<parametersInput.get().size(); pidx++) {
            RealParameter param = parametersInput.get().get(pidx);
            for (int i=0; i<param.getDimension(); i++) {
                if (!indicatorsUsed ||
                        indicatorsInput.get().get(pidx).getValue(i)) {
                    double oldValue = param.getValue(i);
                    double newValue = oldValue*f;
                    if (newValue < param.getLower() || newValue > param.getUpper())
                        return Double.NEGATIVE_INFINITY;

                    param.setValue(i, newValue);
                    logHR += logf;
                }
            }
        }

        // scale parameters inversely
        for (int pidx=0; pidx<parametersInverseInput.get().size(); pidx++) {
            RealParameter param = parametersInverseInput.get().get(pidx);
            for (int i=0; i<param.getDimension(); i++) {
                if (!indicatorsInverseUsed ||
                        indicatorsInverseInput.get().get(pidx).getValue(i)) {
                    double oldValue = param.getValue(i);
                    double newValue = oldValue/f;
                    if (newValue < param.getLower() || newValue > param.getUpper())
                        return Double.NEGATIVE_INFINITY;

                    param.setValue(i, newValue);
                    logHR -= logf;
                }
            }
        }

        // scale tree
        if (!drawEvents) {

            // Note: in this case, events on stem branch can be unrealistic
            // scale node heights and event times
            for (Node node : tree.getNodesAsArray()) {

                node.setHeight(node.getHeight() * f);
                if (!node.isLeaf()) {
                    logHR += logf; // do not count 0 changes!
                }

                for(int i = 0; i < ((AnnotatedNode)node).getEventCount(); i++) {
                    EventNode e = ((AnnotatedNode)node).getEvent(i);
                    double oldTime = e.getHeight();
                    e.setHeight(f * oldTime);
                    if (i > 0) { // do not count the change twice!
                        logHR += logf;
                    }
                }
            }

        } else {

            double newBranchProb;

            // scale node heights only
            for (Node node : tree.getNodesAsArray()) {
                node.setHeight(node.getHeight() * f);
                if (!node.isLeaf()) {
                    logHR += logf; // do not count 0 changes!
                }
            }

            // redraw events (only after new node/ parent heights are set!)
            // if shape or lifetime parameter has changed, the events will be drawn accordingly
            for (int i = 0; i < tree.getNodeCount(); i++) {
                newBranchProb = resampleEvents((AnnotatedNode)tree.getNode(i), true);
                logHR += (oldBranchProbs[i] - newBranchProb);
            }
        }

        // return Hastings ratio
        return logHR;
    }


    // copied from https://github.com/CompEvol/beast2/blob/master/src/beast/base/evolution/operator/ScaleOperator.java
    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
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
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, 1.0 - 1e-8), 1e-8);
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