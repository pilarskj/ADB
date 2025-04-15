package test.adbp;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.util.Randomizer;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class BranchScaleOperatorTest {

    @Test
    public void testSimpleCommands() throws Exception {

        // define tree
        Tree tree = new TreeParser("(((A:10,B:10):6,(C:2,D:2):14):2,(E:7,F:7):11);", true); // nice and easy (quite balanced)
        // Tree tree = new TreeParser("(A:17,(B:14,(C:11,(D:8,(E:5,F:5):3):3):3):3);", true); // hierarchical
        //System.out.println(tree.getInternalNodes());

        double origin = 20;
        double lifetime = 5;
        double branchProportion = 0.8;

        double b = 10;
        double a = lifetime / b;
        //GammaDistribution gammaDist = new GammaDistribution(b, a);

        int branchCount = (int) (branchProportion * tree.getInternalNodeCount());
        int[] internalNodeNumbers = tree.getInternalNodes().stream().mapToInt(Node::getNr).toArray();
        Randomizer.shuffle(internalNodeNumbers);
        int[] randomNodeNumbers = Arrays.copyOfRange(internalNodeNumbers, 0, branchCount);
        System.out.println(Arrays.toString(randomNodeNumbers));

        /*// before manipulation
        for (Node n : tree.getNodesAsArray()) {
            System.out.println("Node number " + n.getNr() + ", height " + n.getHeight());
        } */

        List<Node> reversedNodes = tree.getInternalNodes();
        Collections.reverse(reversedNodes);
        for (Node n : reversedNodes) {
            if (Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getNr())) {
                //double newBranchLength = gammaDist.sample();
                double newBranchLength = lifetime;
                double newHeight;
                if (n.isRoot()) {
                    newHeight = origin - newBranchLength;
                } else {
                    newHeight = n.getParent().getHeight() - newBranchLength;
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

}
