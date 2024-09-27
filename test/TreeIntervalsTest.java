import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.Node;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

public class TreeIntervalsTest {

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
    public void testTreeEdgeTimes() {

        // Define tree
        String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        // Traverse all nodes and compute start/end times of all edges (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);

            // Start time is the node height (or divergence time of the node)
            double startTime = node.getHeight();

            if (!node.isRoot()) {
                // End time is the parent node height (i.e., the divergence time of the parent)
                double endTime = node.getParent().getHeight();

                if (node.isLeaf()) {
                    System.out.println("External edge (Leaf " + node.getID() + "): Start time = " + startTime + ", End time = " + endTime);
                } else {
                    System.out.println("Internal edge (Node " + node.getNr() + "): Start time = " + startTime + ", End time = " + endTime);
                }
            }
        }
    }

    @Test
    public void testBranchingTimes() {

        // Define tree
        String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);


        // Get branching times
        BranchingTimesList B = BranchingTimes.getBranchingTimes(tree);

        System.out.println(Arrays.toString(B.internalStartTimes));
        System.out.println(Arrays.toString(B.internalEndTimes));
        System.out.println(Arrays.toString(B.externalStartTimes));
        System.out.println(Arrays.toString(B.externalEndTimes));
    }
}
