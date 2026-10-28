package adb.operators;

import adb.distribution.ADBTreeDistribution;
import adb.distribution.Parameterization;
import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedTreeParser;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * The conditional event count of AnnotatedTreeOperator (input treeDistribution, single type): the count probabilities,
 * the leaf draw, the agreement of draw and evaluation, and that the proposal is the conditional of the events given
 * the node heights up to the variation of P0 along each branch.
 */
public class ConditionalEventCountTest {

    AnnotatedTreeParser tree;
    RealParameter lifetime, shape, death, sampling;
    Parameterization model;
    ADBTreeDistribution density;

    @BeforeEach
    public void setUp() {
        Randomizer.setSeed(1);
        tree = new AnnotatedTreeParser();
        tree.initByName("newick", AnnotatedTreeOperatorTest.NEWICK);

        lifetime = new RealParameter("1.0");
        shape = new RealParameter("3.0"); // small shape: counts and positions vary most
        death = new RealParameter("0.1");
        sampling = new RealParameter("0.3");
        model = new Parameterization();
        model.initByName("nTypes", 1, "lifetime", lifetime, "shape", shape, "death", death, "sampling", sampling,
                "originTime", 7.859634);
        density = new ADBTreeDistribution();
        density.initByName("tree", tree, "parameterization", model);
        density.calculateLogP();
    }

    EventSampler conditionalSampler() {
        EventSampler operator = new EventSampler();
        operator.initByName("tree", tree, "parameterization", model, "treeDistribution", density,
                "proportionBranches", 1.0, "optimise", false, "weight", 1.0);
        return operator;
    }

    static double sum(double[] logP) {
        double sum = 0;
        for (double x : logP) sum += Math.exp(x);
        return sum;
    }

    static double max(double[] logP) {
        double max = Double.NEGATIVE_INFINITY;
        for (double x : logP) max = Math.max(max, x);
        return max;
    }

    @Test
    public void testCountProbabilitiesSumToOneAndCoverTheTail() {
        for (boolean leaf : new boolean[]{false, true}) {
            for (double length : new double[]{0.05, 0.7, 3.0, 12.0}) {
                for (double k : new double[]{1.0, 3.0, 50.0}) {
                    for (double c : new double[]{0.3, 1.0, 1.7}) {
                        double[] logP = AnnotatedTreeOperator.conditionalCountLogProbabilities(leaf, length, 1.2, k, c);
                        String id = "leaf " + leaf + ", length " + length + ", k " + k + ", c " + c;
                        assertEquals(1.0, sum(logP), 1e-12, id);
                        assertTrue(logP[logP.length - 1] < max(logP) - AnnotatedTreeOperator.COUNT_CUTOFF + 1e-9,
                                "the last term is not negligible: " + id);
                    }
                }
            }
        }
    }

    // with shape 1 and c = 1, the waiting times are exponential, so the number of hidden events n - 1 is Poisson
    // (rate * length), on internal branches and leaves alike
    @Test
    public void testCountIsPoissonForExponentialWaitingTimes() {
        double length = 2.3, lifetime = 0.8;
        PoissonDistribution poisson = new PoissonDistribution(length / lifetime);
        for (boolean leaf : new boolean[]{false, true}) {
            double[] logP = AnnotatedTreeOperator.conditionalCountLogProbabilities(leaf, length, lifetime, 1.0, 1.0);
            for (int n = 1; n <= 15; n++) {
                assertEquals(poisson.logProbability(n - 1), logP[n - 1], 1e-9, "leaf " + leaf + ", n " + n);
            }
        }
    }

    // log Q in log space equals log of commons-math's Q where that does not underflow, and follows the asymptote
    // e^-x x^(a-1) / Gamma(a) (1 + (a - 1) / x) where it does
    @Test
    public void testLogRegularizedGammaQ() {
        for (double a : new double[]{0.5, 1.0, 3.0, 50.0, 100.0}) {
            for (double x = 0.0; x <= 5 * a + 10; x += a / 20 + 0.01) {
                double q = org.apache.commons.math3.special.Gamma.regularizedGammaQ(a, x);
                if (q < 1e-250) continue;
                assertEquals(Math.log(q), AnnotatedTreeOperator.logRegularizedGammaQ(a, x), 1e-9 * Math.max(1, -Math.log(q)),
                        "a " + a + ", x " + x);
            }
            double x = 50 * a + 1000;
            double asymptote = -x + (a - 1) * Math.log(x) - org.apache.commons.math3.special.Gamma.logGamma(a)
                    + Math.log1p((a - 1) / x);
            assertEquals(asymptote, AnnotatedTreeOperator.logRegularizedGammaQ(a, x), 1e-3, "a " + a + ", x " + x);
        }
    }

    // the tip weight of the tilted leaf draw is monotone with value 1 at 0, increasing if (k - 1)(1 - c) > 0 and
    // decreasing otherwise (proof in dev/conditional-event-sampling.md), so max(1, weight(length)) bounds it
    @Test
    public void testTipWeightIsMonotone() {
        for (double k : new double[]{0.5, 1.0, 2.0, 3.0, 10.0, 50.0, 100.0}) {
            for (double c : new double[]{0.05, 0.2, 0.5, 0.9, 1.3, 2.0}) {
                double rate = k / 1.3, tiltedRate = rate * Math.pow(c, 1 / k);
                double direction = Math.signum((k - 1) * (1 - c));
                double previous = 0; // log weight at u = 0
                for (double u = 0.01; u <= 30; u += 0.01) {
                    double logW = AnnotatedTreeOperator.logTipWeight(u, k, rate, tiltedRate);
                    assertTrue(Double.isFinite(logW), "k " + k + ", c " + c + ", u " + u);
                    assertTrue(direction * (logW - previous) >= -1e-9, "not monotone: k " + k + ", c " + c + ", u " + u);
                    if (direction == 0) assertEquals(0.0, logW, 1e-9);
                    previous = logW;
                }
            }
        }
    }

    // the tilted draw gives the exact conditional: the counts follow the count probabilities, the exact density of
    // the sum s integrates to 1 and the drawn s match its mean, also where rejection from the untilted gamma failed
    // (a long leaf with one hidden event, k 3, c 0.2)
    @Test
    public void testLeafWaitingTimesAreExact() {
        double[][] cases = {{10.0, 3.0, 0.2}, {10.0, 100.0, 0.2}, {2.5, 3.0, 1.3}, {3.0, 0.5, 1.6}, {5.0, 50.0, 0.5}};
        double lifetime = 1.0;
        for (double[] x : cases) {
            double length = x[0], k = x[1], c = x[2], scale = lifetime / k;
            String id = "length " + length + ", k " + k + ", c " + c;
            double[] logP = AnnotatedTreeOperator.conditionalCountLogProbabilities(true, length, lifetime, k, c);
            int draws = 20_000, mode = 1;
            for (int n = 1; n <= logP.length; n++) if (logP[n - 1] > logP[mode - 1]) mode = n;
            if (mode == 1) mode = 2; // s needs at least one full waiting time
            int[] counts = new int[logP.length + 1];
            double sum = 0, sumSquares = 0;
            int atMode = 0;
            for (int i = 0; i < draws; i++) {
                double[] w = AnnotatedTreeOperator.drawLeafWaitingTimes(length, lifetime, k, c, logP);
                counts[w.length + 1]++;
                if (w.length + 1 == mode) {
                    double s = 0;
                    for (double v : w) s += v;
                    assertTrue(s > 0 && s < length, id);
                    sum += s;
                    sumSquares += s * s;
                    atMode++;
                }
            }
            for (int n = 1; n <= logP.length; n++) {
                double p = Math.exp(logP[n - 1]);
                double se = Math.sqrt(p * (1 - p) / draws);
                assertEquals(p, (double) counts[n] / draws, 4 * se + 1e-4, id + ", n " + n);
            }

            int points = 200_000;
            double h = length / points, integral = 0, mean = 0;
            for (int i = 0; i < points; i++) {
                double s = (i + 0.5) * h;
                double d = Math.exp(AnnotatedTreeOperator.logLeafSumTarget(s, mode, length, k, scale));
                integral += d * h;
                mean += s * d * h;
            }
            double tolerance = (mode - 1) * k < 1 ? 1e-2 : 1e-4; // the density is infinite at 0 for (n - 1) k < 1
            assertEquals(1.0, integral, tolerance, id);
            if (atMode > 100) {
                double drawnMean = sum / atMode;
                double se = Math.sqrt((sumSquares / atMode - drawnMean * drawnMean) / atMode);
                assertEquals(mean, drawnMean, 4 * se + tolerance, id + ", n " + mode);
            }
        }
    }

    // c = 0 (no extinction risk, so no hidden events): only n = 1
    @Test
    public void testZeroFactorAllowsOnlyOneWaitingTime() {
        for (boolean leaf : new boolean[]{false, true}) {
            double[] logP = AnnotatedTreeOperator.conditionalCountLogProbabilities(leaf, 2.0, 1.0, 3.0, 0.0);
            assertEquals(0.0, logP[0], 1e-12);
            for (int i = 1; i < logP.length; i++) assertEquals(Double.NEGATIVE_INFINITY, logP[i]);
        }
    }

    @Test
    public void testCountSupportBeyondTheCapIsAnError() {
        assertThrows(IllegalStateException.class,
                () -> AnnotatedTreeOperator.conditionalCountLogProbabilities(false, 30_000, 1.0, 3.0, 2.0));
    }

    // the exact mean of the interpolated P0 against a fine midpoint rule, also beyond the grid (origin 7.86)
    @Test
    public void testMeanP0() {
        double[][] intervals = {{0, 0.001}, {0, 1.3}, {0.37, 0.3712}, {2.1, 6.9}, {7.0, 9.0}};
        for (double[] interval : intervals) {
            int points = 100_000;
            double sum = 0;
            for (int i = 0; i < points; i++) {
                sum += density.getP0(0, interval[0] + (i + 0.5) * (interval[1] - interval[0]) / points);
            }
            assertEquals(sum / points, density.getMeanP0(0, interval[0], interval[1]), 1e-8,
                    interval[0] + " to " + interval[1]);
        }
    }

    // the density returned by resampleEvents is the one getBranchProbability evaluates for the same events
    @Test
    public void testDrawDensityMatchesEvaluatedDensity() {
        EventSampler operator = conditionalSampler();
        for (int draw = 0; draw < 50; draw++) {
            for (Node node : tree.getNodesAsArray()) {
                double drawn = operator.resampleEvents((AnnotatedNode) node, true);
                double evaluated = operator.getBranchProbability((AnnotatedNode) node, true);
                assertEquals(drawn, evaluated, 1e-9, "node " + node.getNr());
            }
        }
    }

    // the drawn counts on an internal branch and on a leaf follow the count probabilities
    @Test
    public void testCountFrequencies() {
        EventSampler operator = conditionalSampler();
        for (AnnotatedNode node : new AnnotatedNode[]{(AnnotatedNode) tree.getRoot().getLeft(),
                (AnnotatedNode) tree.getExternalNodes().get(0)}) {
            double[] logP = countLogProbabilities(operator, node);
            int draws = 20_000;
            int[] counts = new int[logP.length + 1];
            for (int i = 0; i < draws; i++) {
                operator.resampleEvents(node, true);
                counts[node.getEventCount()]++;
            }
            for (int n = 1; n <= logP.length; n++) {
                double p = Math.exp(logP[n - 1]);
                double se = Math.sqrt(p * (1 - p) / draws);
                assertEquals(p, (double) counts[n] / draws, 4 * se + 1e-4, "node " + node.getNr() + ", n " + n);
            }
        }
    }

    // the count probabilities of a branch below the root, with c from the operator
    double[] countLogProbabilities(EventSampler operator, AnnotatedNode node) {
        double length = node.getParent().getHeight() - node.getHeight();
        return AnnotatedTreeOperator.conditionalCountLogProbabilities(node.isLeaf(), length, lifetime.getValue(),
                shape.getValue(), operator.hiddenEventFactor(node.getHeight(), node.getParent().getHeight()));
    }

    // log density - log proposal of the tree, minus the variation of log P0 at the hidden events around the branch
    // means, is the same for every redraw: the proposal is the conditional of the events up to that variation
    double residual(EventSampler operator) {
        double residual = density.calculateLogP();
        for (Node n : tree.getNodesAsArray()) {
            AnnotatedNode node = (AnnotatedNode) n;
            residual -= operator.getBranchProbability(node, true);
            double bottom = node.getHeight(), top = node.isRoot() ? model.getOriginTime() : node.getParent().getHeight();
            double mean = density.getMeanP0(0, bottom, top);
            for (int e = 1; e < node.getEventCount(); e++) {
                residual -= Math.log(density.getP0(0, node.getEvent(e).getHeight()) / mean);
            }
        }
        return residual;
    }

    @Test
    public void testProposalIsTheConditionalUpToP0() {
        EventSampler operator = conditionalSampler();
        double before = residual(operator);
        for (int draw = 0; draw < 20; draw++) {
            operator.proposal();
            assertEquals(before, residual(operator), 1e-6, "draw " + draw);
        }
    }

    // solveP0 gives P0 under the current parameter values, while the density keeps its solve until BEAST marks it dirty
    @Test
    public void testSolveP0() {
        double[][] P0 = density.solveP0();
        for (double[] interval : new double[][]{{0, 1.3}, {2.1, 6.9}, {7.0, 9.0}}) {
            assertEquals(density.getMeanP0(0, interval[0], interval[1]), density.getMeanP0(P0, 0, interval[0], interval[1]), 1e-12);
        }
        lifetime.setValue(1.5);
        P0 = density.solveP0();
        ADBTreeDistribution fresh = new ADBTreeDistribution();
        fresh.initByName("tree", tree, "parameterization", model);
        for (double[] interval : new double[][]{{0, 1.3}, {2.1, 6.9}, {7.0, 9.0}}) {
            double mean = fresh.getMeanP0(0, interval[0], interval[1]);
            assertEquals(mean, density.getMeanP0(P0, 0, interval[0], interval[1]), 1e-12);
            assertNotEquals(mean, density.getMeanP0(0, interval[0], interval[1]), 1e-6, "the density keeps the old solve");
        }
    }

    // UpDownEventSampler draws the new events under P0 of the new parameters: its log Hastings ratio is the old
    // densities under P0 of the old parameters, minus the new ones under P0 of a density solved at the new parameters,
    // plus log f for the lifetime
    @Test
    public void testUpDownEventSamplerUsesP0OfTheNewParameters() {
        lifetime.setBounds(0.0, Double.POSITIVE_INFINITY);
        UpDownEventSampler upDown = new UpDownEventSampler();
        upDown.initByName("tree", tree, "parameterization", model, "treeDistribution", density, "up", lifetime,
                "proportionBranches", 1.0, "scaleFactor", 0.3, "optimise", false, "weight", 1.0);
        for (int proposal = 0; proposal < 20; proposal++) {
            density.initAndValidate(); // P0 of the current parameters, as after an evaluation in the MCMC
            double oldLifetime = lifetime.getValue();
            double expected = 0;
            for (Node node : tree.getNodesAsArray()) expected += upDown.getBranchProbability((AnnotatedNode) node, true);

            double logHR = upDown.proposal();

            ADBTreeDistribution fresh = new ADBTreeDistribution();
            fresh.initByName("tree", tree, "parameterization", model);
            EventSampler evaluator = new EventSampler();
            evaluator.initByName("tree", tree, "parameterization", model, "treeDistribution", fresh, "weight", 1.0);
            for (Node node : tree.getNodesAsArray()) expected -= evaluator.getBranchProbability((AnnotatedNode) node, true);
            expected += Math.log(lifetime.getValue() / oldLifetime);
            assertEquals(expected, logHR, 1e-9, "proposal " + proposal);
        }
    }

    // as for UpDownEventSampler, for AnnotatedScale with the lifetime; the Jacobian is log f per internal node and
    // for the lifetime (one type, so no type changes are mapped)
    @Test
    public void testAnnotatedScaleUsesP0OfTheNewParameters() {
        lifetime.setBounds(0.0, Double.POSITIVE_INFINITY);
        AnnotatedScale scale = new AnnotatedScale();
        scale.initByName("tree", tree, "parameterization", model, "treeDistribution", density, "parameter", lifetime,
                "scaleFactor", 0.1, "optimise", false, "weight", 1.0);
        for (int proposal = 0; proposal < 20; proposal++) {
            density.initAndValidate();
            double oldLifetime = lifetime.getValue();
            double expected = 0;
            for (Node node : tree.getNodesAsArray()) expected += scale.getBranchProbability((AnnotatedNode) node, true);

            double logHR = scale.proposal();
            if (logHR == Double.NEGATIVE_INFINITY) continue; // root above the origin, nothing changed

            ADBTreeDistribution fresh = new ADBTreeDistribution();
            fresh.initByName("tree", tree, "parameterization", model);
            EventSampler evaluator = new EventSampler();
            evaluator.initByName("tree", tree, "parameterization", model, "treeDistribution", fresh, "weight", 1.0);
            for (Node node : tree.getNodesAsArray()) expected -= evaluator.getBranchProbability((AnnotatedNode) node, true);
            expected += (tree.getInternalNodeCount() + 1) * Math.log(lifetime.getValue() / oldLifetime);
            assertEquals(expected, logHR, 1e-9, "proposal " + proposal);
        }
    }
}
