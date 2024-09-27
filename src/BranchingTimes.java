import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.Arrays;

// Object: stores branching times
class BranchingTimesList {
    double[] internalStartTimes;
    double[] internalEndTimes;
    double[] externalStartTimes;
    double[] externalEndTimes;

    public BranchingTimesList(double[] int_s, double[] int_e, double[] ext_s, double[] ext_e) {
        this.internalStartTimes = int_s;
        this.internalEndTimes = int_e;
        this.externalStartTimes = ext_s;
        this.externalEndTimes = ext_e;
    }
}


public class BranchingTimes {

    // Function for getting branching times from a Tree
    public static BranchingTimesList getBranchingTimes(Tree tree) {

        // create empty arrays to store start and end times of branches
        double[] int_s = new double[0];
        double[] int_e = new double[0];
        double[] ext_s = new double[0];
        double[] ext_e = new double[0];

        // traverse all nodes and compute start/end times of all branches (backwards in time)
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);

            // start time is the node height (or divergence time of the node)
            double startTime = node.getHeight();

            if (!node.isRoot()) {
                // end time is the parent node height (i.e., the divergence time of the parent)
                double endTime = node.getParent().getHeight();

                // if node is a tip, add to external branches, otherwise, add to internal branches
                if (!node.isLeaf()) {
                    // resize the arrays and accommodate the new elements
                    int_s = Arrays.copyOf(int_s, int_s.length + 1);
                    int_e = Arrays.copyOf(int_e, int_e.length + 1);
                    int_s[int_s.length - 1] = startTime;
                    int_e[int_e.length - 1] = endTime;
                } else {
                    ext_s = Arrays.copyOf(ext_s, ext_s.length + 1);
                    ext_e = Arrays.copyOf(ext_e, ext_e.length + 1);
                    ext_s[ext_s.length - 1] = startTime;
                    ext_e[ext_e.length - 1] = endTime;
                }
            }
        }

        BranchingTimesList B = new BranchingTimesList(int_s, int_e, ext_s, ext_e);
        return B;
    }
}
