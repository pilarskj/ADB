import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.Node;

import org.junit.jupiter.api.Test;

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
        String newick = "((D:5.0,C:4.0):6.0,(A:1.0,B:2.0):3.0):0.0;";

        // Parse tree
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        // Traverse all nodes and compute start/end times of all edges
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);

            // End time is the node height (or divergence time of the node)
            double endTime = node.getHeight();

            if (!node.isRoot()) {
                // Start time is the parent node height (i.e., the divergence time of the parent)
                double startTime = node.getParent().getHeight();

                if (node.isLeaf()) {
                    System.out.println("External edge (Leaf " + node.getID() + "): Start time = " + startTime + ", End time = " + endTime);
                } else {
                    System.out.println("Internal edge (Node " + node.getNr() + "): Start time = " + startTime + ", End time = " + endTime);
                }
            }
        }
    }
}
