package adb.operators;

import adb.distribution.ADBTreeDistribution;
import adb.distribution.Parameterization;
import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedNode.EventNode;
import adb.tree.AnnotatedTree;
import adb.util.Utils;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math3.distribution.GammaDistribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.ContinuedFraction;

import static org.apache.commons.math3.special.Gamma.regularizedGammaP;
import static org.apache.commons.math3.special.Gamma.regularizedGammaQ;


/**
 * Base class for operators that change an AnnotatedTree together with its hidden events.
 * Operators keep type changes (see scaleEvents) and redraw hidden events between them with resampleEvents;
 * getBranchProbability evaluates the same draw density for current events, and the difference enters the log Hastings ratio.
 * Per type segment, the hidden events are drawn as follows:
 * - default: the count from a discretised normal (renewal CLT), the positions from Dirichlet(shape); a leaf at time 0
 *   by a backward renewal draw;
 * - with treeDistribution (single type): the count from its conditional given the branch length, with P0 replaced
 *   by its branch mean, and the positions from their exact conditional (on leaves by drawLeafWaitingTimes). P0 is that
 *   of the current state, so an operator that changes parameters before it redraws must set proposedP0.
 */
@Description("Base class of operators that propose changes to an annotated tree and its hidden events")
public abstract class AnnotatedTreeOperator extends TreeOperator {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);

    public Input<ADBTreeDistribution> treeDistributionInput = new Input<>("treeDistribution",
            "ADB tree distribution (single type): if given, the event counts are drawn from their conditional "
                    + "given the branch length instead of the normal approximation");

    protected Parameterization parameterization;
    protected ADBTreeDistribution treeDistribution;
    // P0 of the conditional count under parameters changed within the proposal (from treeDistribution.solveP0);
    // null: P0 of the current state
    protected double[][] proposedP0;

    // cut-off and cap of the count support
    static final double COUNT_CUTOFF = 40;
    static final int MAX_COUNT = 100_000;


    @Override
    public void initAndValidate() {

        Tree tree = treeInput.get();
        if (!(tree instanceof AnnotatedTree)) {
            throw new IllegalArgumentException("Attempted to initialise annotated tree operator with regular tree.");
        }

        parameterization = parameterizationInput.get();
        treeDistribution = treeDistributionInput.get();
        if (treeDistribution != null && parameterization.getNTypes() != 1) {
            throw new IllegalArgumentException("The conditional event count (treeDistribution) is only implemented for a single type");
        }
    }


    /**
     * Whether the model allows the division of the parent into the two daughters, given the type at the parent
     * and the type at the bottom of each daughter's branch. Moves that create an impossible division have zero
     * likelihood, so rejecting them early only saves the cost of resampling and of the likelihood.
     */
    protected boolean isDivisionPossible(Node parent, Node left, Node right) {
        return parameterization.getDivisionFactor(AnnotatedNode.getTypeOf(parent),
                ((AnnotatedNode) left).getInitialType(), ((AnnotatedNode) right).getInitialType()) > 0.0;
    }

    
    /**
     * Redraws the hidden events above the node between its type changes and returns the log density of the draw.
     * The top is the parent, or the origin for the root (without origin nothing is drawn). If drawEventCount is false,
     * the counts are kept (except on a leaf at time 0, which uses the default draw).
     */
    protected double resampleEvents(AnnotatedNode node, boolean drawEventCount) {

        double logp = 0;

        double tMax;
        if (node.isRoot()) {
            if (parameterization.getOriginTime() == null) {
                node.getEvent(0).setHeight(node.getHeight()); // keep the node event in sync with the node height
                return 0.0; // no stem: the root has no hidden events, so there is nothing to draw
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

        // return probability
        return logp;
    }


    // Draws the hidden events of one type segment [start, end) into events and returns the log density of the draw.
    private double sampleEventsPerType(AnnotatedNode node, double start, double end, int type, int nEvents, boolean drawEventCount, ArrayList<EventNode> events) {

        double logp = 0;

        double lifetime = parameterization.getLifetime(type);
        double shape = parameterization.getShape(type);
        double scale = lifetime / shape;

        boolean conditional = treeDistribution != null && drawEventCount;

        if (node.isLeaf() && start == 0.0 && conditional) {
            double c = hiddenEventFactor(start, end);
            double[] logP = conditionalCountLogProbabilities(true, end, lifetime, shape, c);
            double[] waitingTimes = drawLeafWaitingTimes(end, lifetime, shape, c, logP); // from the top down
            int n = waitingTimes.length + 1;
            logp += logP[n - 1];
            if (n > 1) {
                double t = end;
                for (double w : waitingTimes) t -= w;
                logp += logLeafSumTarget(end - t, n, end, shape, scale) + logDirichlet(waitingTimes, shape);
                for (int i = waitingTimes.length - 1; i > 0; i--) { // the events in ascending order
                    events.add(new EventNode(type, t));
                    t += waitingTimes[i];
                }
                events.add(new EventNode(type, t));
            }

        } else if (node.isLeaf() && start == 0.0) {
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
            if (conditional) {
                double[] logP = conditionalCountLogProbabilities(false, length, lifetime, shape, hiddenEventFactor(start, end));
                eventCount = drawCount(logP);
                logp += logP[eventCount - 1];
            } else if (drawEventCount) {
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


    /**
     * Log density with which resampleEvents would draw the current hidden events above the node
     * (the same segments and the same distributions, but evaluated instead of sampled).
     */
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


    // Log density of the hidden events of one type segment, starting at event index ix, mirrors sampleEventsPerType
    private double getProbabilityPerType(AnnotatedNode node, int ix, double start, double end, int type, int eventCount, boolean drawEventCount) {

        double logp = 0;

        double lifetime = parameterization.getLifetime(type);
        double shape = parameterization.getShape(type);
        double scale = lifetime / shape;

        double waitingTime;
        boolean conditional = treeDistribution != null && drawEventCount;

        if (node.isLeaf() && start == 0.0 && conditional) {
            double[] logP = conditionalCountLogProbabilities(true, end, lifetime, shape, hiddenEventFactor(start, end));
            if (eventCount > logP.length) return Double.NEGATIVE_INFINITY;
            logp += logP[eventCount - 1];
            if (eventCount > 1) {
                double s = end - node.getEvent(1).getHeight();
                logp += logLeafSumTarget(s, eventCount, end, shape, scale) + logDirichlet(node, 1, eventCount - 1, end, shape);
            }

        } else if (node.isLeaf() && start == 0.0) {
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
            if (conditional) {
                double[] logP = conditionalCountLogProbabilities(false, length, lifetime, shape, hiddenEventFactor(start, end));
                if (eventCount > logP.length) return Double.NEGATIVE_INFINITY;
                logp += logP[eventCount - 1];
            } else if (drawEventCount) {
                double mean = length / lifetime;
                double variance = ((shape * (scale * scale)) * length) / Math.pow(lifetime, 3);
                double sd = Math.sqrt(variance);
                logp += Utils.logProbabilityDiscretizedNormal(eventCount, mean, sd);
            }
            logp += logDirichlet(node, ix, eventCount, end, shape);
        }

        return logp;
    }


    // Dirichlet log density of the eventCount waiting times from event ix to end
    private static double logDirichlet(AnnotatedNode node, int ix, int eventCount, double end, double shape) {
        double[] waitingTimes = new double[eventCount];
        for (int i = 0; i < eventCount; i++) {
            double upper = i + 1 < eventCount ? node.getEvent(ix + i + 1).getHeight() : end;
            waitingTimes[i] = upper - node.getEvent(ix + i).getHeight();
        }
        return logDirichlet(waitingTimes, shape);
    }

    // log density of waiting times whose proportions of their sum are Dirichlet(shape), as drawn by
    // Utils.distributeDirichletProbability
    private static double logDirichlet(double[] waitingTimes, double shape) {
        int n = waitingTimes.length;
        double length = 0;
        for (double w : waitingTimes) length += w;
        double logDirichlet = 0;
        for (double w : waitingTimes) logDirichlet += Math.log(w / length);
        return -(n - 1) * Math.log(length) + Gamma.logGamma(n * shape) - n * Gamma.logGamma(shape) + (shape - 1) * logDirichlet;
    }


    // factor c = (1 - d)(2 sym + asym) P0 of a hidden event, with the mean P0 on [start, end]
    double hiddenEventFactor(double start, double end) {
        return (1 - parameterization.getDeath(0))
                * (2 * parameterization.getSymTransition(0, 0) + parameterization.getAsymTransition(0, 0))
                * (proposedP0 == null ? treeDistribution.getMeanP0(0, start, end)
                        : treeDistribution.getMeanP0(proposedP0, 0, start, end));
    }

    /**
     * Normalised log probabilities of the number n = 1, 2, ... of waiting times (index n - 1) on a branch,
     * for constant P0 (c from hiddenEventFactor): p(n) ∝ c^(n-1) Gamma(n shape)(length);
     * for a leaf at height 0, p(n) ∝ c^(n-1) P(S_{n-1} < length < S_n).
     * The support depends only on the arguments (terms down to e^-COUNT_CUTOFF of the largest, at least
     * 2 max(1, c) length / lifetime + 5), so draw and evaluation agree.
     */
    static double[] conditionalCountLogProbabilities(boolean leaf, double length, double lifetime, double shape, double c) {
        double rate = shape / lifetime;
        int minCount = (int) Math.ceil(2 * Math.max(1, c) * length / lifetime) + 5;
        double[] terms = new double[Math.min(minCount, MAX_COUNT)];
        double max = Double.NEGATIVE_INFINITY;
        int size = 0;
        while (true) {
            int n = size + 1;
            if (n > MAX_COUNT) {
                throw new IllegalStateException("The count probabilities need more than " + MAX_COUNT
                        + " terms on a branch of length " + length + " (lifetime " + lifetime + ", c " + c + ")");
            }
            // c^0 = 1 also for c = 0, where only n = 1 is possible
            double term = (n == 1 ? 0 : (n - 1) * Math.log(c))
                    + (leaf ? logLeafCountTerm(n, length, shape, rate) : logGammaDensity(length, n * shape, rate));
            if (size == terms.length) terms = Arrays.copyOf(terms, 2 * size);
            terms[size++] = term;
            max = Math.max(max, term);
            if (n >= minCount && term < max - COUNT_CUTOFF) break;
        }
        double sum = 0;
        for (int i = 0; i < size; i++) sum += Math.exp(terms[i] - max);
        double logNorm = max + Math.log(sum);
        double[] logP = new double[size];
        for (int i = 0; i < size; i++) logP[i] = terms[i] - logNorm;
        return logP;
    }

    // log P(S_{n-1} < length < S_n), S_m the sum of m waiting times
    static double logLeafCountTerm(int n, double length, double shape, double rate) {
        double x = rate * length;
        if (n == 1) return logRegularizedGammaQ(shape, x);
        double below = regularizedGammaP((n - 1) * shape, x);
        // upper tails are more accurate when the lower ones are close to 1
        double z = below > 0.5 ? regularizedGammaQ(n * shape, x) - regularizedGammaQ((n - 1) * shape, x)
                : below - regularizedGammaP(n * shape, x);
        return z > 0 ? Math.log(z) : Double.NEGATIVE_INFINITY;
    }

    private static double logGammaDensity(double x, double shape, double rate) {
        return shape * Math.log(rate) + (shape - 1) * Math.log(x) - rate * x - Gamma.logGamma(shape);
    }

    // draws n from conditionalCountLogProbabilities
    private static int drawCount(double[] logP) {
        double u = Randomizer.nextDouble();
        double cumulative = 0;
        for (int i = 0; i < logP.length; i++) {
            cumulative += Math.exp(logP[i]);
            if (u < cumulative) return i + 1;
        }
        return logP.length; // u above the rounded sum
    }

    // exact log density of the sum s of the n - 1 full waiting times above a leaf, given n:
    // Gamma((n - 1) shape)(s) (1 - F(length - s)) / P(S_{n-1} < length < S_n)
    static double logLeafSumTarget(double s, int n, double length, double shape, double scale) {
        double rate = 1 / scale;
        return logGammaDensity(s, (n - 1) * shape, rate) + logRegularizedGammaQ(shape, rate * (length - s))
                - logLeafCountTerm(n, length, shape, rate);
    }

    /**
     * Draws the full waiting times above a leaf (from the top down) exactly from their conditional
     * ∝ c^(n-1) ∏ g(w_i) (1 - F(u)), with u the tip interval, restricted to the counts with positive probability in logP.
     * Exponential tilting gives c g(w) = g'(w) e^(-theta w), with g' the gamma density with rate' = rate c^(1/shape) and
     * theta = rate - rate'. So a backward renewal draw with g' (the waiting times until the first that does not fit)
     * differs from the target only by tipWeight(u) = e^(theta u) (1 - F(u)) / (1 - F'(u)), and is accepted with
     * probability tipWeight(u) / max(1, tipWeight(length)). The bound holds because tipWeight is monotone with
     * tipWeight(0) = 1.
     */
    static double[] drawLeafWaitingTimes(double length, double lifetime, double shape, double c, double[] logP) {
        if (c == 0) return new double[0]; // no hidden events
        double rate = shape / lifetime;
        double tiltedRate = rate * Math.pow(c, 1 / shape);
        double logBound = Math.max(0, logTipWeight(length, shape, rate, tiltedRate));
        double[] waitingTimes = new double[8];
        while (true) {
            int count = 0;
            double sum = 0;
            while (true) {
                double w = Randomizer.nextGamma(shape, tiltedRate);
                if (sum + w >= length) break;
                if (count == waitingTimes.length) waitingTimes = Arrays.copyOf(waitingTimes, 2 * count);
                waitingTimes[count++] = w;
                sum += w;
            }
            if (count + 1 > logP.length || logP[count] == Double.NEGATIVE_INFINITY) continue; // outside the support
            if (Math.log(Randomizer.nextDouble()) < logTipWeight(length - sum, shape, rate, tiltedRate) - logBound) {
                return Arrays.copyOf(waitingTimes, count);
            }
        }
    }

    // log tipWeight(u) = theta u + log(1 - F(u)) - log(1 - F'(u)), with F, F' the gamma CDFs with rate and tiltedRate
    static double logTipWeight(double u, double shape, double rate, double tiltedRate) {
        return (rate - tiltedRate) * u + logRegularizedGammaQ(shape, rate * u) - logRegularizedGammaQ(shape, tiltedRate * u);
    }

    // log Q(a, x) without underflow: log(1 - P) below a + 1, otherwise the continued fraction of
    // Gamma.regularizedGammaQ, Q = e^(-x) x^a / Gamma(a) / cf(x), in log space
    static double logRegularizedGammaQ(double a, double x) {
        if (x < a + 1) return Math.log1p(-regularizedGammaP(a, x));
        ContinuedFraction cf = new ContinuedFraction() {
            @Override
            protected double getA(int n, double x) { return 2.0 * n + 1.0 - a + x; }

            @Override
            protected double getB(int n, double x) { return n * (a - n); }
        };
        return -x + a * Math.log(x) - Gamma.logGamma(a) - Math.log(cf.evaluate(x, 1e-15, Integer.MAX_VALUE));
    }


    /**
     * Maps the events above node i from the branch [i.height, oldTop] affinely onto [i.height, newTop],
     * keeping the node itself fixed. Used by the topology operators when the parent of i changes height.
     * Only the type changes need to be kept (the other hidden events are resampled afterwards), so the return value
     * is the log Jacobian of this deterministic map for the m type changes: m * log(scale).
     * It is 0 for a single type, and the branch is then fully described by resampleEvents.
     */
    protected static double scaleEvents(AnnotatedNode i, double oldTop, double newTop) {
        return scaleEvents(i, i.getHeight(), oldTop, i.getHeight(), newTop, false);
    }

    /**
     * Maps the events above node i affinely from the branch [oldBase, oldTop] onto [newBase, newTop], where newBase
     * is the new height of i, so the order of the events and their branch are kept. Returns the log Jacobian
     * m * log(scale) of this deterministic map, where m counts the type changes, or all events above the node if
     * allEvents is true (when none of them is resampled afterwards).
     */
    protected static double scaleEvents(AnnotatedNode i, double oldBase, double oldTop, double newBase, double newTop,
                                        boolean allEvents) {
        double scale = (newTop - newBase) / (oldTop - oldBase);
        i.getEvent(0).setHeight(newBase);

        int mapped = 0;
        for (int e = 1; e < i.getEventCount(); e++) {
            EventNode event = i.getEvent(e);
            event.setHeight(newBase + (event.getHeight() - oldBase) * scale);
            if (allEvents || event.getType() != i.getEvent(e - 1).getType()) {
                mapped++;
            }
        }
        return mapped == 0 ? 0.0 : mapped * Math.log(scale);
    }

    /**
     * Multiplies entry index of every parameter by factor, or every entry if index is -1. Returns the share of the
     * log Hastings ratio of a scale move with a kernel symmetric on the log scale, log(factor) per scaled entry,
     * or minus infinity if a new value is outside the bounds of its parameter.
     */
    protected static double scaleParameters(List<RealParameter> parameters, int index, double factor) {
        int count = 0;
        for (RealParameter param : parameters) {
            int from = index < 0 ? 0 : index;
            int to = index < 0 ? param.getDimension() : index + 1;
            for (int i = from; i < to; i++) {
                double newValue = param.getValue(i) * factor;
                if (newValue < param.getLower() || newValue > param.getUpper()) return Double.NEGATIVE_INFINITY;
                param.setValue(i, newValue);
                count++;
            }
        }
        return count * Math.log(factor);
    }

}
