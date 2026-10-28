package adb.operators;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import adb.tree.AnnotatedNode;

/**
 * Adapted from the Uniform operator of BEAST. With an origin, the root is moved too, between its children and the
 * origin; its third branch is the stem. The type changes on the three branches at the node are mapped affinely
 * onto their new branches, and the hidden events of the three branches are redrawn (with treeDistribution from their
 * conditional). The height draw is symmetric, so the log Hastings ratio is the log Jacobian of the type changes
 * plus the old minus the new event density of the three branches.
 */
@Description("Moves an internal node height uniformly between its children and parent (the origin for the root, which moves only with an origin), and redraws the hidden events of its three branches")
public class AnnotatedUniform extends AnnotatedTreeOperator {

    @Override
    public double proposal() {
        final Tree tree = (Tree) InputUtil.get(treeInput, this);
        final int nodeCount = tree.getNodeCount();

        // the root moves only below an origin
        final Double origin = parameterization.getOriginTime();

        // abort if no movable internal nodes
        if (tree.getInternalNodeCount() == 1 && origin == null)
            return Double.NEGATIVE_INFINITY;

        // randomly select internal node
        Node node;
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while ((node.isRoot() && origin == null) || node.isLeaf());

        // draw new height
        final double upper = node.isRoot() ? origin : node.getParent().getHeight();
        final double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;

        // calculate pre-move branch density
        AnnotatedNode[] branches = {(AnnotatedNode) node, (AnnotatedNode) node.getLeft(), (AnnotatedNode) node.getRight()};
        double logHR = 0;
        for (AnnotatedNode branch : branches) {
            logHR += getBranchProbability(branch, true);
        }

        // set new height
        final double oldValue = node.getHeight();
        node.setHeight(newValue);

        // scale type changes and draw new events along branches
        logHR += scaleEvents(branches[0], oldValue, upper, newValue, upper, false);
        logHR += scaleEvents(branches[1], oldValue, newValue) + scaleEvents(branches[2], oldValue, newValue);
        for (AnnotatedNode branch : branches) {
            logHR -= resampleEvents(branch, true);
        }

        return logHR;
    }

}
