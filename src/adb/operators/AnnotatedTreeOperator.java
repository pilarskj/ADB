package adb.operators;

import adb.distribution.Parameterization;
import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedTree;
import adb.util.Utils;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.Arrays;

@Description("This operator generates proposals for an AnnotatedTree.")
public abstract class AnnotatedTreeOperator extends TreeOperator {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.OPTIONAL);

    protected Parameterization parameterization;


    @Override
    public void initAndValidate() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);
        if (!(tree instanceof AnnotatedTree)) {
            throw new IllegalArgumentException("Attempted to initialise annotated tree operator with regular tree.");
        }

        parameterization = parameterizationInput.get();
    }


    // TODO: extend to multi-type branches
    // only makes sense for internal branches
    protected double redistributeEvents(AnnotatedNode node) {

        double shape = parameterization.getShape(0);
        double origin;

        if (node.isLeaf()) {
            return Double.NEGATIVE_INFINITY;
        } else {
            if (node.isRoot()) {
                if (parameterization.getOriginTime() == null) {
                    return 0.0; // TODO: correct?
                } else {
                    origin = parameterization.getOriginTime();
                }
            } else {
                origin = node.getParent().getHeight();
            }
            return redistributeEvents(node, origin, shape);
        }
    }


    private double redistributeEvents(AnnotatedNode node, double origin, double shape) {
        // TODO: align with distributeEvents function in AnnotatedNode?
        int n = node.getEventCount();
        double[] segmentLengths = new double[n];

        for (int i = 0; i < (n - 1); i++) {
            segmentLengths[i] = node.getEvent(i + 1).getHeight() - node.getEvent(i).getHeight();
        }
        segmentLengths[n - 1] = origin - node.getEvent(n - 1).getHeight();
        double logp = Utils.distributeDirichletProbability(segmentLengths, shape);

        double height = node.getHeight();
        for (int i = 0; i < (n - 1); i++) {
            height = height + segmentLengths[i];
            node.getEvent(i + 1).setHeight(height);
        }

        return logp;
    }


    protected double resampleEvents(AnnotatedNode node) {

        double lifetime = parameterization.getLifetime(0);
        double shape = parameterization.getShape(0);
        double rate = shape / lifetime;
        GammaDistribution gammaDist = new GammaDistribution(shape, 1/rate);

        // clear current events
        int n = node.getEventCount();
        node.clearEvents();

        double tMin = node.getHeight();
        double tMax;
        double t;
        double logp = 0;

        if (node.isRoot()) {
            if (parameterization.getOriginTime() == null) {
                return 0.0; // TODO: correct?
            } else {
                tMax = parameterization.getOriginTime();
            }
        } else {
            tMax = node.getParent().getHeight();
        }

        if (node.isLeaf()) {
            t = tMax;
            double waitingTime;
            while (t > tMin) { // tMin should be 0 in this case!
                // determine the next hidden event
                waitingTime = Randomizer.nextGamma(shape, rate);
                t -= waitingTime;
                if (t > tMin) {
                    // add new hidden event
                    node.addEvent(new AnnotatedNode.EventNode(0, t), 1);
                    logp += gammaDist.logDensity(waitingTime);
                }
            }
            // final time interval
            double remainder;
            if (node.getEventCount() > 1) {
                remainder = node.getEvent(1).getHeight();
            } else {
                remainder = tMax;
            }
            logp += Math.log(1 - gammaDist.cumulativeProbability(remainder));

        } else {
            // sample events
            double length = tMax - tMin;
            //int eventCount = n;
            int eventCount = (int) Randomizer.nextPoisson(length/lifetime);
            logp += (new PoissonDistribution(length/lifetime)).logProbability(eventCount);
            if (eventCount > 1) {
                double[] segmentLengths = new double[eventCount];
                Arrays.fill(segmentLengths, length/eventCount);
                // distribute
                logp += Utils.distributeDirichletProbability(segmentLengths, shape);
                t = tMin;
                for (int i = 1; i < eventCount; i++) {
                    t += segmentLengths[i];
                    node.addEvent(new AnnotatedNode.EventNode(0, t), true);
                }
            }
        }

        // return probability
        return logp;
    }


    protected double getBranchProbability(AnnotatedNode node) {
        double lifetime = parameterization.getLifetime(0);
        double shape = parameterization.getShape(0);
        double scale = lifetime / shape;
        GammaDistribution gammaDist = new GammaDistribution(shape, scale);

        double logp = 0;
        double tMax;
        if (node.isRoot()) {
            if (parameterization.getOriginTime() == null) {
                return 0.0;
            } else {
                tMax = parameterization.getOriginTime();
            }
        } else {
            tMax = node.getParent().getHeight();
        }

        int eventCount = node.getEventCount();
        double waitingTime;

        if (node.isLeaf()) {
            if (eventCount == 1) {
                waitingTime = tMax;
            } else {
                for (int i = 2; i < eventCount; i++) {
                    waitingTime = node.getEvent(i).getHeight() - node.getEvent(i - 1).getHeight();
                    logp += gammaDist.logDensity(waitingTime);
                }
                waitingTime = tMax - node.getEvent(eventCount - 1).getHeight();
                logp += gammaDist.logDensity(waitingTime);
                // remainder
                waitingTime = node.getEvent(1).getHeight();
            }
            logp += Math.log(1 - gammaDist.cumulativeProbability(waitingTime));

        } else {
            double length = tMax - node.getHeight();
            logp += (new PoissonDistribution(length/lifetime)).logProbability(eventCount);
            double logDirichlet = 0;
            for (int i = 1; i < eventCount; i++) {
                waitingTime = node.getEvent(i).getHeight() - node.getEvent(i - 1).getHeight();
                logDirichlet += Math.log(waitingTime/length);
            }
            waitingTime = tMax - node.getEvent(eventCount - 1).getHeight();
            logDirichlet += Math.log(waitingTime/length);
            logp += -(eventCount - 1) * Math.log(length) +
                    Gamma.logGamma(eventCount * shape) - eventCount * Gamma.logGamma(shape) + (shape - 1) * logDirichlet;
        }

        return logp;
    }

}
