package test.adb;

import adb.Branch;
import adb.BranchList;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.Node;

import org.junit.jupiter.api.Test;

public class BranchTest {

    @Test
    public void testTreeClass() {

        // define tree
        String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";
        double origin = 15;

        // parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("newick", newick, "IsLabelledNewick", true, "adjustTipHeights", false);

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

        // define tree
        String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";
        double origin = 15;

        // parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("newick", newick, "IsLabelledNewick", true, "adjustTipHeights", false);

        BranchList branches = new BranchList(tree, origin);

        for (Branch branch : branches.listBranches()) {
            System.out.println("Branch from node " + branch.startNode + " to node " + branch.endNode);
            System.out.println("Start time: " + branch.startTime + ", End time: " + branch.endTime);
            System.out.println("Branch type: " + branch.branchType);
            System.out.println("Branch index: " + branch.branchIndex + ", left: " + branch.leftIndex + ", right: " + branch.rightIndex);
            System.out.println("Mode: " + branch.branchMode);
            System.out.println();
        }
    }


    @Test
    public void testTreeTraversalColoured() {

        // define tree
        // from R simulator
        // branch-typed:
        String newick = "((((1[&type=0]:1.378074977,(11[&type=1]:0.9265429956)[&type=1]:0.4515319812)[&type=1]:0.5929298741)[&type=0]:2.088429726,(((3[&type=0]:0.4513672607)[&type=1]:1.262419189,((6[&type=0]:0.05632331623)[&type=1]:0.7255498933)[&type=1]:0.9319132401)[&type=1]:0.2056271944)[&type=0]:2.140020933)[&type=0]:2.073580344)[&type=0]:1.866985079;";
        // tip-typed:
        //String newick = "((1[&type=0]:1.378074977,11[&type=1]:1.378074977):2.6813596,(3[&type=0]:1.71378645,6[&type=0]:1.71378645):2.345648127):3.940565423;";
        double originTime = 8;
        int originType = 0;

        // parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("newick", newick, "IsLabelledNewick", true, "adjustTipHeights", true); // due to rounding errors

        /* // traverse all nodes (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) {
                System.out.println("Node Nr " + node.getNr() + ", ID " + node.getID() + ", Type " + node.getMetaData("type") +
                        ", Height " + node.getHeight() + ", Left " + node.getLeft() + ", Right " + node.getRight());
            } else {
                System.out.println("Node Nr " + node.getNr() + ", ID " + node.getID() + ", Type " + node.getMetaData("type") +
                        ", Height " + node.getHeight() + ", Left " + node.getLeft().getNr() + ", Right " + node.getRight().getNr());
            }
        } */

        // check BranchList
        BranchList branches = new BranchList(tree, originTime);

        for (Branch branch : branches.listBranches()) {
            System.out.println();
            System.out.println("Branch from node " + branch.startNode + " to node " + branch.endNode);
            System.out.println("Start time: " + branch.startTime + ", End time: " + branch.endTime);
            System.out.println("Branch type: " + branch.branchType);
            System.out.println("Branch index: " + branch.branchIndex + ", left: " + branch.leftIndex + ", right: " + branch.rightIndex);
            System.out.println("Mode: " + branch.branchMode);
        }
    }
}