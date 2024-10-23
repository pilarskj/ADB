package test.adbp;

import adbp.Branch;
import adbp.BranchList;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.Node;

import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;


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
                System.out.println(i + " External edge (Leaf " + node.getID() + "): Start time = " + startTime + ", End time = " + endTime);
            } else {
                System.out.println(i + " Internal edge (Node " + node.getNr() + "): Start time = " + startTime + ", End time = " + endTime);
            }
        }
    }

    @Test
    public void testTreeTraversal() {

        // Define tree
        //String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";
        String newick = "(C:11.0,(A:8.0,B:8.0):3.0):0.0;";

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        double t_or = 12;

        BranchList branches = new BranchList();
        branches.traverseTree(tree, t_or);
        branches.assignBranchIndices();
        branches.getChildBranches();

        for (Branch branch : branches.getBranches()) {
            System.out.println("Branch from node " + branch.startNode + " to node " + branch.endNode);
            System.out.println("Start time: " + branch.startTime + ", End time: " + branch.endTime);
            System.out.println("Edge index: " + branch.branchIndex);
            System.out.println("Left downstream edge index: " + branch.leftIndex);
            System.out.println("Right downstream edge index: " + branch.rightIndex);
            System.out.println("Type: " + branch.branchType);
            System.out.println();
        }
        /*
        // traverse all nodes and compute start/end times of all branches (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            System.out.println("Node Nr " + node.getNr() + ", ID " + node.getID() + ", Height " + node.getHeight()
                    + ", Left " + node.getLeft() + ", Right " + node.getRight());
        }
        */
    }
}
