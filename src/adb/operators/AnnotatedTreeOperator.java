package adb.operators;

import adb.distribution.Parameterization;
import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedNode.EventNode;
import adb.tree.AnnotatedTree;
import adb.util.Utils;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math3.distribution.GammaDistribution;

import java.util.ArrayList;
import java.util.Arrays;


@Description("This operator generates proposals for an AnnotatedTree.")
public abstract class AnnotatedTreeOperator extends TreeOperator {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.OPTIONAL);

    protected Parameterization parameterization;


    @Override
    public void initAndValidate() {

        Tree tree = treeInput.get();
        if (!(tree instanceof AnnotatedTree)) {
            throw new IllegalArgumentException("Attempted to initialise annotated tree operator with regular tree.");
        }

        parameterization = parameterizationInput.get();
    }


    // TODO: add option to operate on one type only?
    protected double resampleEvents(AnnotatedNode node, boolean drawEventCount) {

        double logp = 0;

        double tMax;
        if (node.isRoot()) {
            if (parameterization.getOriginTime() == null) {
                return 0.0; // TODO: correct?
            } else {
                tMax = parameterization.getOriginTime();
            }
        } else {
            tMax = node.getParent().getHeight();
        }

        int nTotal = node.getEventCount();
        double start = node.getHeight();
        int type = node.getType();
        ArrayList<EventNode> events = new ArrayList<>();
        events.add(new EventNode(type, start));

        int ix = 1;
        int nEvents = 1;
        double end;
        while (ix < nTotal) {
            if (node.getEvent(ix).getType() == type) { // no type change
                nEvents++;
            } else {
                end = node.getEvent(ix).getHeight();
                logp += sampleEventsPerType(node, start, end, type, nEvents, drawEventCount, events);
                // reset counters
                type = node.getEvent(ix).getType();
                start = node.getEvent(ix).getHeight();
                nEvents = 1;
            }
            ix++;
        }
        logp += sampleEventsPerType(node, start, tMax, type, nEvents, drawEventCount, events);

        node.setEvents(events);
        node.makeAllDirty(Tree.IS_DIRTY); // TODO: necessary here?

        // return probability
        return logp;
    }


    // Draw hidden events along a (single-type) branch segment according to lifetime distribution
    private double sampleEventsPerType(AnnotatedNode node, double start, double end, int type, int nEvents, boolean drawEventCount, ArrayList<EventNode> events) {

        double logp = 0;

        double lifetime = parameterization.getLifetime(type);
        double shape = parameterization.getShape(type);
        double scale = lifetime / shape;

        if (node.isLeaf() && start == 0.0) {
            GammaDistribution gammaDist = new GammaDistribution(null, shape, scale);
            double t = end;
            double waitingTime;
            while (t > start) {
                // determine the next hidden event
                waitingTime = Randomizer.nextGamma(shape, 1/scale);
                if (t - waitingTime > start) {
                    t -= waitingTime;
                    // add new hidden event
                    events.add(1, new EventNode(type, t));
                    logp += gammaDist.logDensity(waitingTime);
                } else {
                    // calculate probability of remaining time
                    logp += Math.log(1 - gammaDist.cumulativeProbability(t));
                    break;
                }
            }

        } else {
            if (start > node.getHeight()) {
                events.add(new EventNode(type, start)); // copy initial event
            }

            double length = end - start;
            int eventCount;
            if (drawEventCount) {
                // use CLT for renewal processes: approximate number of events with discretized normal distribution
                double mean = length / lifetime;
                double variance = ((shape * (scale * scale)) * length) / Math.pow(lifetime, 3);
                double sd = Math.sqrt(variance);
                eventCount = Utils.sampleDiscretizedNormal(mean, sd);
                logp += Utils.logProbabilityDiscretizedNormal(eventCount, mean, sd);
            } else {
                eventCount = nEvents;
            }

            if (eventCount > 1) {
                double[] segmentLengths = new double[eventCount];
                Arrays.fill(segmentLengths, length/eventCount);
                // distribute
                logp += Utils.distributeDirichletProbability(segmentLengths, shape);
                double t = start;
                for (int i = 0; i < eventCount - 1; i++) { // last segment corresponds to time until next branching event
                    t += segmentLengths[i];
                    events.add(new EventNode(type, t));
                }
            }
        }

        return logp;
    }


    protected double getBranchProbability(AnnotatedNode node, boolean drawEventCount) {

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

        int nTotal = node.getEventCount();
        double start = node.getHeight();
        int type = node.getType();

        int ix = 1;
        int nEvents = 1;
        double end;
        int ixTypeChange = 0;
        while (ix < nTotal) {
            if (node.getEvent(ix).getType() == type) { // no type change
                nEvents++;
            } else {
                end = node.getEvent(ix).getHeight();
                logp += getProbabilityPerType(node, ixTypeChange, start, end, type, nEvents, drawEventCount);
                // reset counters
                type = node.getEvent(ix).getType();
                start = node.getEvent(ix).getHeight();
                nEvents = 1;
                ixTypeChange = ix;
            }
            ix++;
        }
        logp += getProbabilityPerType(node, ixTypeChange, start, tMax, type, nEvents, drawEventCount);

        // return probability
        return logp;
    }


    private double getProbabilityPerType(AnnotatedNode node, int ix, double start, double end, int type, int eventCount, boolean drawEventCount) {

        double logp = 0;

        double lifetime = parameterization.getLifetime(type);
        double shape = parameterization.getShape(type);
        double scale = lifetime / shape;

        double waitingTime;

        if (node.isLeaf() && start == 0.0) {
            GammaDistribution gammaDist = new GammaDistribution(null, shape, scale);
            if (eventCount == 1) {
                waitingTime = end;
            } else {
                for (int i = 2; i < eventCount; i++) {
                    waitingTime = node.getEvent(i).getHeight() - node.getEvent(i - 1).getHeight();
                    logp += gammaDist.logDensity(waitingTime);
                }
                waitingTime = end - node.getEvent(eventCount - 1).getHeight();
                logp += gammaDist.logDensity(waitingTime);
                // remainder
                waitingTime = node.getEvent(1).getHeight();
            }
            logp += Math.log(1 - gammaDist.cumulativeProbability(waitingTime));

        } else {
            double length = end - start;
            if (drawEventCount) {
                double mean = length / lifetime;
                double variance = ((shape * (scale * scale)) * length) / Math.pow(lifetime, 3);
                double sd = Math.sqrt(variance);
                logp += Utils.logProbabilityDiscretizedNormal(eventCount, mean, sd);
            }
            double logDirichlet = 0;
            for (int i = ix + 1; i < ix + eventCount; i++) {
                waitingTime = node.getEvent(i).getHeight() - node.getEvent(i - 1).getHeight();
                logDirichlet += Math.log(waitingTime / length);
            }
            waitingTime = end - node.getEvent(ix + eventCount - 1).getHeight();
            logDirichlet += Math.log(waitingTime / length);
            logp += -(eventCount - 1) * Math.log(length) +
                    Gamma.logGamma(eventCount * shape) - eventCount * Gamma.logGamma(shape) + (shape - 1) * logDirichlet;
        }

        return logp;
    }

}
