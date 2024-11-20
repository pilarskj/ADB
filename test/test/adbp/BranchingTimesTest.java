package test.adbp;

import adbp.Branch;
import adbp.BranchList;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.Node;

import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

public class BranchingTimesTest {

    @Test
    public void testTreeIntervals() {

        // Define tree
        String newick = "((D:5.0,C:4.0):6.0,(A:1.0,B:2.0):3.0):0.0;";

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        // Instantiate TreeIntervals for this tree
        TreeIntervals treeIntervals = new TreeIntervals(tree);

        // Calculate intervals and output information
        System.out.println("Tree Intervals:");

        // Iterate over the intervals and print them
        for (int i = 0; i < treeIntervals.getIntervalCount(); i++) {
            double interval = treeIntervals.getInterval(i);
            System.out.println("Interval " + i + ": " + interval);
        } // gives waiting times between any event to occur in the tree (not branch lengths!)
    }


    @Test
    public void testBranchingTimes() {
        /*
        // Define tree
        String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);
        */

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double t_or = 15;

        // traverse all nodes and compute start/end times of all branches (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);

            // start time is the node height (or divergence time of the node)
            double startTime = node.getHeight();

            double endTime;
            if (node.isRoot()) {
                // end time is the tree origin
                endTime = t_or;
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
        // String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";
        String newick = "(C:11.0,(A:8.0,B:8.0):3.0):0.0;";
        double origin = 12;

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        BranchList branches = new BranchList(tree, origin, false);
        // branches.assignRandomTypes(2);

        for (Branch branch : branches.listBranches()) {
        // Branch branch = branches.getBranchByIndex(3);
            System.out.println("Branch from node " + branch.startNode + " to node " + branch.endNode);
            System.out.println("Start time: " + branch.startTime + ", End time: " + branch.endTime);
            System.out.println("Branch index: " + branch.branchIndex + ", left: " + branch.leftIndex + ", right: " + branch.rightIndex);
            System.out.println("Mode: " + branch.branchMode);
            // System.out.println("Node type: " + branch.nodeType);
            System.out.println();
        }

        //System.out.println(branches.countExternalBranches());
        //System.out.println(branches.countInternalBranches());
    }

    @Test
    public void testTreeTraversalColoured() {

        // Define tree
        // from R simulator
        // String newick = "((4[&type=1]:9.195045362,(2[&type=1]:0.6232736482,3[&type=1]:0.6232736482)[&type=1]:8.571771714)[&type=1]:0.4465947482,(1[&type=1]:9.505634401,(6[&type=1]:0.3789614042,7[&type=0]:0.3789614042)[&type=1]:9.126672996)[&type=1]:0.1360057095)[&type=0];";
        // from Java simulator
        String newick = "((2[&type=1]:6.287294724563209,(4[&type=0]:1.6450230674511555,5[&type=0]:1.6450230674511555)[&type=0]:4.642271657112053)[&type=0]:0.24321442405762905,1[&type=1]:6.530509148620838)[&type=1]:0.0;";
        double origin = 10;

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", true); // due to rounding errors

        /*
        // simulated tree with 26 tips (typed)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree_typed.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);
        double origin = 15;
        */

        // traverse all nodes (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) {
                System.out.println("Node Nr " + node.getNr() + ", ID " + node.getID() + ", Type " + node.getMetaData("type") +
                        ", Height " + node.getHeight() + ", Left " + node.getLeft() + ", Right " + node.getRight());
            } else {
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
