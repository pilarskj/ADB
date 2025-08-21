package test.adb;

import adb.operators.SampleBranchesOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

public class SampleBranchesOperatorTest {

    @Test
    public void testSimpleCommands() throws Exception {

        // define tree
        Tree tree = new TreeParser("(((A:10,B:10):6,(C:2,D:2):14):2,(E:7,F:7):11);", true); // nice and easy (quite balanced)
        // Tree tree = new TreeParser("(A:17,(B:14,(C:11,(D:8,(E:5,F:5):3):3):3):3);", true); // hierarchical
        //System.out.println(tree.getInternalNodes());

        /*// before manipulation
        for (Node n : tree.getNodesAsArray()) {
            System.out.println("Node number " + n.getNr() + ", height " + n.getHeight());
        }*/

        double origin = 20;
        double lifetime = 5;
        double branchProportion = 0.5;

        double b = 100;
        double a = lifetime / b;
        GammaDistribution gammaDist = new GammaDistribution(b, a);

        int branchCount = (int) (branchProportion * tree.getInternalNodeCount());
        List<Node> internalNodes = tree.getInternalNodes();
        int[] internalNodeNumbers = tree.getInternalNodes().stream().mapToInt(Node::getNr).toArray();
        Randomizer.shuffle(internalNodeNumbers);
        int[] randomNodeNumbers = Arrays.copyOfRange(internalNodeNumbers, 0, branchCount);
        System.out.println(Arrays.toString(randomNodeNumbers));

        /*List<Node> reversedNodes = tree.getInternalNodes();
        Collections.reverse(reversedNodes);
        for (Node n : reversedNodes) { */
        for (int i = internalNodes.size(); i > 0; i--) {
            Node n = internalNodes.get(i - 1);
            if (Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getNr())) {
                double newBranchLength = gammaDist.sample();
                double newHeight;
                if (n.isRoot()) {
                    newHeight = origin - newBranchLength;
                } else {
                    newHeight = n.getParent().getHeight() - newBranchLength;
                }
                // check whether new branch lengths fit in the process
                if (newHeight < 0) {
                    System.out.println("order violated");
                }
                if (!Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getLeft().getNr()) && (newHeight < n.getLeft().getHeight())) {
                    System.out.println("order violated");
                }
                if (!Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getRight().getNr()) && (newHeight < n.getRight().getHeight())) {
                    System.out.println("order violated");
                }
                n.setHeight(newHeight);
            }
        }

        // after manipulation
        for (Node n : tree.getNodesAsArray()) {
            System.out.println("Node number " + n.getNr() + ", height " + n.getHeight());
        }
        System.out.println(tree);
    }

    @Test
    public void testSampleBranches() throws Exception {

        // define tree
        Tree tree = new TreeParser("(((A:10,B:10):6,(C:2,D:2):14):2,(E:7,F:7):11);", true); // nice and easy (quite balanced)
        // Tree tree = new TreeParser("(A:17,(B:14,(C:11,(D:8,(E:5,F:5):3):3):3):3);", true); // hierarchical
        //System.out.println(tree.getInternalNodes());
        System.out.println(tree);

        SampleBranchesOperator op = new SampleBranchesOperator();
        op.setInputValue("tree", tree);
        op.setInputValue("lifetime", new RealParameter("5"));
        op.setInputValue("shapeInteger", new IntegerParameter("10"));
        op.setInputValue("origin", new RealParameter("20"));

        double hr = op.proposal();
        System.out.println(hr);
        System.out.println(tree);
    }

    @Test
    public void testSubtrees() throws Exception {
        Tree tree = new TreeParser("(((A:10,B:10):6,(C:2,D:2):14):2,(E:7,F:7):11);", true);
        Node n = tree.getNode(8);
        List<Node> subtree = n.getAllChildNodesAndSelf();
        for (Node m : subtree) {
            System.out.println("Node number " + m.getNr() + ", height " + m.getHeight());
        }
    }

}
