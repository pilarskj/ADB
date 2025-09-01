package test.adb;

import adb.Branch;
import adb.BranchList;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.Node;

import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

public class BranchTest {

    @Test
    public void testBranchingTimes() {

        // Define tree
        String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("newick", newick, "IsLabelledNewick", true, "adjustTipHeights", false);

        double origin = 15;

        // traverse all nodes and compute start/end times of all branches (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);

            // start time is the node height (or divergence time of the node)
            double startTime = node.getHeight();

            double endTime;
            if (node.isRoot()) {
                // end time is the tree origin
                endTime = origin;
            } else {
                // end time is the parent node height (i.e., the divergence time of the parent)
                endTime = node.getParent().getHeight();
            }

            if (node.isLeaf()) {
                System.out.println(i + " External branch (Leaf " + node.getID() + "): Start time = " + startTime + ", End time = " + endTime);
            } else {
                System.out.println(i + " Internal branch (Node " + node.getNr() + "): Start time = " + startTime + ", End time = " + endTime);
            }
        }
    }

    @Test
    public void testTreeTraversal() {

        // Define tree
        String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";
        double origin = 15;

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("newick", newick, "IsLabelledNewick", true, "adjustTipHeights", false);

        BranchList branches = new BranchList(tree, origin, false);
        // branches.assignRandomTypes(2);

        for (Branch branch : branches.listBranches()) {
            // Branch branch = branches.getBranchByIndex(3);
            System.out.println("Branch from node " + branch.startNode + " to node " + branch.endNode);
            System.out.println("Start time: " + branch.startTime + ", End time: " + branch.endTime);
            System.out.println("Branch index: " + branch.branchIndex + ", left: " + branch.leftIndex + ", right: " + branch.rightIndex);
            System.out.println("Mode: " + branch.branchMode);
            System.out.println("Node type: " + branch.nodeType);
            System.out.println();
        }

        // System.out.println(branches.countExternalBranches());
        // System.out.println(branches.countInternalBranches());
    }

    @Test
    public void testTreeTraversalColoured() {

        // Define tree
        // from R simulator
        String newick = "((4[&type=1]:9.195045362,(2[&type=1]:0.6232736482,3[&type=1]:0.6232736482)[&type=1]:8.571771714)[&type=1]:0.4465947482,(1[&type=1]:9.505634401,(6[&type=1]:0.3789614042,7[&type=0]:0.3789614042)[&type=1]:9.126672996)[&type=1]:0.1360057095)[&type=0];";
        // String newick = "(((14[&type=0]:4.072008045,(10[&type=1]:2.499314806,(8[&type=1]:0.7573405951,9[&type=1]:0.7573405951)[&type=0]:1.741974211)[&type=0]:1.572693238)[&type=0]:1.809062642,(18[&type=0]:3.741049754,(1[&type=0]:1.558300887,2[&type=0]:1.558300887)[&type=0]:2.182748867)[&type=0]:2.140020933)[&type=0]:2.251944234,(3[&type=0]:6.059434577,(22[&type=1]:3.990901207,((26[&type=1]:0.3738130627,27[&type=1]:0.3738130627)[&type=1]:1.501512376,(20[&type=1]:0.8660709291,21[&type=1]:0.8660709291)[&type=1]:1.009254509)[&type=0]:2.115575768)[&type=0]:2.06853337)[&type=0]:2.073580344):1.866985079;";
        // from Java simulator
        // String newick = "((2[&type=1]:6.287294724563209,(4[&type=0]:1.6450230674511555,5[&type=0]:1.6450230674511555)[&type=0]:4.642271657112053)[&type=0]:0.24321442405762905,1[&type=1]:6.530509148620838)[&type=1]:0.0;";
        double origin = 10;
        double origin_type = 0;

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("newick", newick, "IsLabelledNewick", true, "adjustTipHeights", true); // due to rounding errors

        // traverse all nodes (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) {
                System.out.println("Node Nr " + node.getNr() + ", ID " + node.getID() + ", Type " + node.getMetaData("type") +
                        ", Height " + node.getHeight() + ", Left " + node.getLeft() + ", Right " + node.getRight());
            } else {
                if (node.isRoot() & node.getMetaData("type") == null) {
                    node.setMetaData("type", origin_type);
                }
                System.out.println("Node Nr " + node.getNr() + ", ID " + node.getID() + ", Type " + node.getMetaData("type") +
                        ", Height " + node.getHeight() + ", Left " + node.getLeft().getNr() + ", Right " + node.getRight().getNr());
            }
        }

        // check BranchList
        BranchList branches = new BranchList(tree, origin, true);

        for (Branch branch : branches.listBranches()) {
            System.out.println();
            System.out.println("Branch from node " + branch.startNode + " to node " + branch.endNode);
            System.out.println("Start time: " + branch.startTime + ", End time: " + branch.endTime);
            System.out.println("Branch index: " + branch.branchIndex + ", left: " + branch.leftIndex + ", right: " + branch.rightIndex);
            System.out.println("Mode: " + branch.branchMode);
            System.out.println("Node type: " + branch.nodeType);
        }
    }
}