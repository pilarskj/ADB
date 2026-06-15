package adb.operators;

import adb.distribution.Parameterization;
import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedNode.EventNode;
import adb.tree.AnnotatedTree;
import adb.util.Utils;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import org.apache.commons.math3.distribution.GammaDistribution;

import java.util.Arrays;


public class EventSampler extends TreeOperator {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);

    protected Parameterization parameterization;

    @Override
    public void initAndValidate() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);
        if (!(tree instanceof AnnotatedTree)) {
            throw new IllegalArgumentException("Attempted to initialise annotated tree operator with regular tree.");
        }

        parameterization = parameterizationInput.get();

    }

    @Override
    public double proposal() {
        return 0;
    }


    // TODO: extend to multi-type branches
    protected double resampleEvents(AnnotatedNode node) {

        double lifetime = parameterization.getLifetime(0);
        double shape = parameterization.getShape(0);
        double rate = shape / lifetime;
        GammaDistribution gammaDist = new GammaDistribution(shape, 1/rate);

        // clear current events
        node.clearEvents();

        double tMin = node.getHeight();
        double tMax;
        double t;
        double prob = 0;

        if (node.isLeaf()) {
            tMax = node.getParent().getHeight();
            t = tMax;
            double duration;
            while (t > tMin) {
                // determine the next hidden event
                duration = Randomizer.nextGamma(shape, rate);
                t -= duration;
                if (t > tMin) {
                    // add new hidden event
                    node.addEvent( new EventNode(0, t), 1);
                    prob += gammaDist.logDensity(duration);
                }
            }
            // final time interval
            double remainder = node.getEvent(1).getHeight();
            prob += Math.log(1 - gammaDist.cumulativeProbability(remainder));

        } else {
            if (node.isRoot()) {
                if (parameterization.getOriginTime() == null) {
                    return 1; // TODO: correct?
                } else {
                    tMax = parameterization.getOriginTime();
                }
            } else {
                tMax = node.getParent().getHeight();
            }
            double length = tMax - tMin;
            int eventCount = (int)(length / lifetime);
            double[] segmentLengths = new double[eventCount + 1];
            Arrays.fill(segmentLengths, length / (eventCount + 1));
            // distribute // TODO: align with distributeEvents function?
            double[] newSegmentLengths = Utils.distributeDirichlet(segmentLengths, shape);
            t = tMin;
            double duration;
            for (int i = 0; i <= eventCount; i++) {
                duration = newSegmentLengths[i];
                t += duration;
                node.addEvent( new EventNode(0, t), true);
                prob += gammaDist.logDensity(duration);
            }
        }

        // return branch probability
        return prob;

    }


    protected double getBranchProb(AnnotatedNode node) {
        double lifetime = parameterization.getLifetime(0);
        double shape = parameterization.getShape(0);
        double scale = lifetime / shape;
        GammaDistribution gammaDist = new GammaDistribution(shape, scale);

        double prob = 0;
        double tMax;
        if (node.isRoot()) {
            if (parameterization.getOriginTime() == null) {
                return 1;
            } else {
                tMax = parameterization.getOriginTime();
            }
        } else {
            tMax = node.getParent().getHeight();
        }

        double tUpper = tMax;
        double tLower;
        int i = node.getEventCount() - 1;
        while (i > 0) {
            tLower = node.getEvents().get(i).getHeight();
            prob += gammaDist.logDensity(tUpper - tLower);
            tUpper = tLower;
            i--;
        }

        if (node.isLeaf()) {
            prob += Math.log(1 - gammaDist.cumulativeProbability(tUpper));
        } else {
            prob += gammaDist.logDensity(tUpper - node.getHeight());
        }

        return prob;
    }

}
