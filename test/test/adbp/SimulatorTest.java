package test.adbp;

import adbp.Simulator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.Test;

public class SimulatorTest {

    @Test
    public void testTreeSimulator() {

        double origin = 10;
        int originType = 0;
        double[] a = {2, 3};
        double[] b = {1, 2};
        double[] d = {0.2, 0.1};
        double rho = 0.5;
        double[][] Xsi_as = {{0, 0.1},
                {0.2, 0}};
        double[][] Xsi_s = {{0.5, 0.3},
                {0.1, 0.5}};

        Simulator sim = new Simulator();
        //Tree tree = sim.simulateCompleteTree(origin, originType, a, b, d, Xsi_as, Xsi_s);
        //tree = sim.pruneTree(tree, rho);
        Tree tree = sim.simulateTree(origin, originType, a, b, d, Xsi_as, Xsi_s, rho);

        // traverse all nodes
        if (tree != null) {
            for (int i = 0; i < tree.getNodeCount(); i++) {
                Node node = tree.getNode(i);
                if (node.isLeaf()) {
                    System.out.println("Node Nr " + node.getNr() + ", Type " + node.getMetaData("type") +
                            ", Height " + node.getHeight() + ", Left " + node.getLeft() + ", Right " + node.getRight() +
                            ", isRoot " + node.isRoot());
                } else {
                    System.out.println("Node Nr " + node.getNr() + ", Type " + node.getMetaData("type") +
                            ", Height " + node.getHeight() + ", Left " + node.getLeft().getNr() + ", Right " + node.getRight().getNr() +
                            ", isRoot " + node.isRoot());
                }
            }
            // convert to Newick format
            System.out.println(tree.getRoot().toNewick());
        }
    }


    @Test
    public void testTreePruning() {

        //String newick = "((A:6,B:6):4,(C:8,(D:5,E:5):3):2):0;";
        String newick = "(((A:6,B:6):4,(C:8,(D:5,E:5):3):2):2,(F:4,(G:6,H:2):4):2):0;";

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        // traverse all nodes and compute start/end times of all branches (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) {
                System.out.println("Node Nr " + node.getNr() + ", Type " + node.getMetaData("type") +
                        ", Height " + node.getHeight() + ", Left " + node.getLeft() + ", Right " + node.getRight());
            } else {
                System.out.println("Node Nr " + node.getNr() + ", Type " + node.getMetaData("type") +
                        ", Height " + node.getHeight() + ", Left " + node.getLeft().getNr() + ", Right " + node.getRight().getNr());
            }
        }

        System.out.println("\n");
        Simulator sim = new Simulator();
        Tree prunedTree = sim.pruneTree(tree, 0.3);

        if (prunedTree != null) {
            // traverse all nodes and compute start/end times of all branches (backwards in time)
            for (int i = 0; i < prunedTree.getNodeCount(); i++) {
                Node node = prunedTree.getNode(i);
                if (node.isLeaf()) {
                    System.out.println("Node Nr " + node.getNr() + ", Type " + node.getMetaData("type") +
                            ", Height " + node.getHeight() + ", Left " + node.getLeft() + ", Right " + node.getRight() +
                            ", isRoot " + node.isRoot());
                } else {
                    System.out.println("Node Nr " + node.getNr() + ", Type " + node.getMetaData("type") +
                            ", Height " + node.getHeight() + ", Left " + node.getLeft().getNr() + ", Right " + node.getRight().getNr() +
                            ", isRoot " + node.isRoot());
                }
            }
            // convert to Newick format
            System.out.println(tree.getRoot().toNewick());
        }
    }

}
