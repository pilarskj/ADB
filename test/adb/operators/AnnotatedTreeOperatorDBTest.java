package adb.operators;

import adb.distribution.ADBTreeDistribution;
import adb.distribution.Parameterization;
import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedTree;
import adb.tree.AnnotatedTreeInitialiser;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.speciation.YuleModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.DirectSimulator;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Collections;
import java.util.List;
import java.util.Random;

// adapted from https://github.com/NicoNeureiter/targetedbeast/blob/main/test/targetedbeast/WilsonBaldingDBTest.java

/**
 * Detailed-balance test for the targeted tree operators.
 *
 * <p>Specialises {@link DetailedBalanceTest} for {@link Tree} states sampled from
 * a {@link YuleModel} prior, with tree-level quantisation functions
 * (root child height, tree length, imbalance). Each {@code @Test} method builds a
 * {@link TreeTrial} with one operator/edge-weight combination, defined inline.
 */
public class AnnotatedTreeOperatorDBTest extends DetailedBalanceTest<Tree> {

    protected static final int NUM_SAMPLES = 20_000;
    protected static final int NUM_TAXA = 10;
    protected static final double BIRTH_DIFF_RATE = 1.0;
    protected static Parameterization parameterization;

    @BeforeClass
    public static void setUpClass() {
        parameterization = getParameterization();
    }

    @Override
    protected int getNumSamples() {
        return NUM_SAMPLES;
    }

    @Override
    protected List<StateMapper<Tree>> getStateMappers() {
        return List.of(
            new StateMapper<Tree>("RootFirstChildHeight", tree -> {
                List<Node> rootChildren = tree.getRoot().getChildren();
                double firstChildHeight = Math.max(rootChildren.get(0).getHeight(), rootChildren.get(1).getHeight());
                return String.format("%.0f", 2 * firstChildHeight);
            }),
            new StateMapper<Tree>("TreeLength", tree -> {
                double treeLength = TreeUtils.getTreeLength(tree, tree.getRoot());
                return String.format("%.0f", 2 * treeLength);
            }),
            new StateMapper<Tree>("TreeImbalance", tree -> String.valueOf(rootImbalance(tree.getRoot())))
        );
    }

    @Test
    public void testYule() throws Exception {
        testDetailedBalance(new TreeTrial() {
            @Override
            protected Distribution getTargetPrior(Tree tree, EventDensitySupport eventDensity) {
                return new YuleWithHiddenEventDensity(getYulePrior(tree), eventDensity);
            }

            @Override
            protected TreeOperator getOperator(Tree tree) {
                return getAnnotatedScale(tree);  // change here operator for testing
            }
        });
    }

    @Ignore("ADBTreeDistribution cannot directly sample independent annotated trees for this harness yet.")
    @Test
    public void testADBTreeDistribution() throws Exception {
        testDetailedBalance(new TreeTrial() {
            @Override
            protected Distribution getTargetPrior(Tree tree, EventDensitySupport eventDensity) {
                ADBTreeDistribution prior = new ADBTreeDistribution();
                prior.initByName("tree", tree, "parameterization", parameterization);
                return prior;
            }

            @Override
            protected TreeOperator getOperator(Tree tree) {
                return getAnnotatedScale(tree); // change here operator for testing
            }
        });
    }

    /**
     * A {@link Trial} over annotated trees: a plain Yule-distributed branching
     * tree is repeatedly redrawn and converted to a fresh annotated tree before
     * the operator is applied.
     */
    abstract class TreeTrial implements Trial<Tree> {
        final Tree branchingTree = new Tree();
        final AnnotatedTree annotatedTree = new AnnotatedTreeInitialiser();
        final TreeDistribution samplingPrior;
        final EventDensitySupport eventDensity;
        final Distribution targetPrior;
        final DirectSimulator simulator = new DirectSimulator();
        final TreeOperator operator;

        TreeTrial() {
            branchingTree.initByName("taxonset", getTaxonSet(NUM_TAXA));
            samplingPrior = getYulePrior(branchingTree);
            simulator.initByName("distribution", samplingPrior, "nSamples", 1);
            annotatedTree.initByName("branchingTree", branchingTree, "parameterization", parameterization);
            eventDensity = new EventDensitySupport();
            eventDensity.initByName("tree", annotatedTree, "parameterization", parameterization, "weight", 1.0);
            eventDensity.resampleAllEvents(annotatedTree);
            targetPrior = getTargetPrior(annotatedTree, eventDensity);
            operator = getOperator(annotatedTree);
        }

        protected abstract Distribution getTargetPrior(Tree tree, EventDensitySupport eventDensity);

        protected abstract TreeOperator getOperator(Tree tree);

        @Override
        public AnnotatedTree nextSample() throws Exception {
            simulator.run();
            AnnotatedTreeInitialiser sample = new AnnotatedTreeInitialiser();
            sample.initByName("branchingTree", branchingTree, "parameterization", parameterization);
            annotatedTree.assignFrom(sample);
            eventDensity.resampleAllEvents(annotatedTree);
            return annotatedTree;
        }

        @Override
        public Distribution prior() {
            return targetPrior;
        }

        @Override
        public double proposal() {
            return operator.proposal();
        }
    }

    static class YuleWithHiddenEventDensity extends Distribution {
        final TreeDistribution yulePrior;
        final EventDensitySupport eventDensity;

        YuleWithHiddenEventDensity(TreeDistribution yulePrior, EventDensitySupport eventDensity) {
            this.yulePrior = yulePrior;
            this.eventDensity = eventDensity;
        }

        @Override
        public double calculateLogP() {
            logP = yulePrior.calculateLogP() + eventDensity.calculateLogP();
            return logP;
        }

        @Override
        public List<String> getArguments() {
            return Collections.emptyList();
        }

        @Override
        public List<String> getConditions() {
            return Collections.emptyList();
        }

        @Override
        public void sample(State state, Random random) {
            throw new UnsupportedOperationException();
        }
    }

    static class EventDensitySupport extends AnnotatedTreeOperator {
        double calculateLogP() {
            Tree tree = treeInput.get();
            double logP = 0.0;
            for (Node node : tree.getNodesAsArray()) {
                logP += getBranchProbability((AnnotatedNode) node, true);
            }
            return logP;
        }

        void resampleAllEvents(AnnotatedTree tree) {
            for (Node node : tree.getNodesAsArray()) {
                resampleEvents((AnnotatedNode) node, true);
            }
        }

        @Override
        public double proposal() {
            throw new UnsupportedOperationException();
        }
    }

    protected TreeOperator getAnnotatedWilsonBalding(Tree tree) {
        AnnotatedWilsonBalding operator = new AnnotatedWilsonBalding();
        operator.initByName("tree", tree, "parameterization", parameterization, "weight", 1.0);
        return operator;
    }

    protected TreeOperator getAnnotatedNarrowExchange(Tree tree) {
        AnnotatedExchange operator = new AnnotatedExchange();
        operator.initByName("tree", tree, "parameterization", parameterization, "isNarrow", true, "weight", 1.0);
        return operator;
    }

    protected TreeOperator getAnnotatedWideExchange(Tree tree) {
        AnnotatedExchange operator = new AnnotatedExchange();
        operator.initByName("tree", tree, "parameterization", parameterization, "isNarrow", false, "weight", 1.0);
        return operator;
    }

    protected TreeOperator getAnnotatedSubtreeSlide(Tree tree) {
        AnnotatedSubtreeSlide operator = new AnnotatedSubtreeSlide();
        operator.initByName("tree", tree, "parameterization", parameterization, "size", 10, "weight", 1.0);
        return operator;
    }

    protected TreeOperator getAnnotatedScale(Tree tree) {
        AnnotatedScale operator = new AnnotatedScale();
        operator.initByName("tree", tree, "parameterization", parameterization, "drawEvents", true, "weight", 1.0);
        return operator;
    }

    protected TaxonSet getTaxonSet(int numTaxa) {
        TaxonSet taxonSet = new TaxonSet();
        for (int i = 0; i < numTaxa; i++) {
            taxonSet.initByName("taxon", new Taxon(String.valueOf(i)));
        }
        return taxonSet;
    }

    protected TreeDistribution getYulePrior(Tree tree) {
        YuleModel treePrior = new YuleModel();
        treePrior.initByName("tree", tree, "birthDiffRate", String.valueOf(BIRTH_DIFF_RATE));
        return treePrior;
    }

    static public Parameterization getParameterization() {
        Parameterization parameterization = new Parameterization();
        parameterization.initByName("nTypes", 1,
                "lifetime", new RealParameter("0.1"),
                "shape", new IntegerParameter("10"),
                "death", new RealParameter("0.1"),
                "sampling", new RealParameter("0.1")
                //"originTime", 10.0
        );
        return parameterization;
    }
}
