package adb.operators;

import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedNode.EventNode;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

// modified version of beast.base.evolution.operator.Exchange
public class AnnotatedExchange extends AnnotatedTreeOperator {

    public Input<Boolean> isNarrowInput = new Input<>("isNarrow", "if true (default) a narrow exchange is performed, otherwise a wide exchange", true);

    @Override
    public double proposal() {

        final Tree tree = (Tree) InputUtil.get(treeInput, this);

        double logHR;

        if (isNarrowInput.get()) {
            logHR = narrow(tree);
        } else {
            logHR = wide(tree);
        }

        return logHR;
    }

    private int isg(final Node n) {
        return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    private int sisg(final Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }


    public double narrow(final Tree tree) {

        boolean exchangeEvents = true;
        double iOldBranchProb = 0, uncleOldBranchProb = 0, iNewBranchProb, uncleNewBranchProb, logHR = 0;

        final int internalNodes = tree.getInternalNodeCount();
        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        Node grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        while (grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf()) {
            grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        }

        Node parentIndex = grandParent.getLeft();
        Node uncle = grandParent.getRight();
        if (parentIndex.getHeight() < uncle.getHeight()) {
            parentIndex = grandParent.getRight();
            uncle = grandParent.getLeft();
        }

        int validGP = 0;
        {
            for (int i = internalNodes + 1; i < 1 + 2*internalNodes; ++i) {
                validGP += isg(tree.getNode(i));
            }
        }

        int c2 = sisg(parentIndex) + sisg(uncle);

        Node i = (Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight());

        // TODO: extend resampling events to multi-type version
        if (parameterization != null && parameterization.getNTypes() == 1) {
            iOldBranchProb = getBranchProbability((AnnotatedNode)i, true);
            uncleOldBranchProb = getBranchProbability((AnnotatedNode)uncle, true);
            exchangeEvents = false;
        }

        exchangeNodes(i, uncle, parentIndex, grandParent, exchangeEvents);

        if (parameterization != null && parameterization.getNTypes() == 1) {
            // resample events along branch and calculate probability
            iNewBranchProb = resampleEvents((AnnotatedNode)i, true);
            logHR += (iOldBranchProb - iNewBranchProb);
            uncleNewBranchProb = resampleEvents((AnnotatedNode)uncle, true);
            logHR += (uncleOldBranchProb - uncleNewBranchProb);
        }

        int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);
        logHR += Math.log((float)validGP/validGPafter);

        return logHR;
    }


    public double wide(final Tree tree) {

        boolean exchangeEvents = true;
        double iOldBranchProb = 0, jOldBranchProb = 0, iNewBranchProb, jNewBranchProb, logHR = 0;

        final int nodeCount = tree.getNodeCount();
        Node i = tree.getRoot();
        while (i.isRoot()) {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
        }

        Node j = i;
        while (j.getNr() == i.getNr() || j.isRoot()) {
            j = tree.getNode(Randomizer.nextInt(nodeCount));
        }

        final Node p = i.getParent();
        final Node jP = j.getParent();

        if ((p != jP) && (i != jP) && (j != p)
                && (j.getHeight() < p.getHeight())
                && (i.getHeight() < jP.getHeight())) {

            // TODO: extend resampling events to multi-type version
            if (parameterization != null && parameterization.getNTypes() == 1) {
                iOldBranchProb = getBranchProbability((AnnotatedNode)i, true);
                jOldBranchProb = getBranchProbability((AnnotatedNode)j, true);
                exchangeEvents = false;
            }

            exchangeNodes(i, j, p, jP, exchangeEvents);

            // all the nodes on the path from i/j to the common ancestor of i/j parents had a topology change,
            // so they need to be marked FILTHY.
            if (markCladesInput.get()) {
                Node iup = p;
                Node jup = jP;
                while (iup != jup) {
                    if( iup.getHeight() < jup.getHeight() ) {
                        assert !iup.isRoot();
                        iup = iup.getParent();
                        iup.makeDirty(Tree.IS_FILTHY);
                    } else {
                        assert !jup.isRoot();
                        jup = jup.getParent();
                        jup.makeDirty(Tree.IS_FILTHY);
                    }
                }
            }

            if (parameterization != null && parameterization.getNTypes() == 1) {
                // resample events along branch and calculate probability
                iNewBranchProb = resampleEvents((AnnotatedNode)i, true);
                logHR += (iOldBranchProb - iNewBranchProb);
                jNewBranchProb = resampleEvents((AnnotatedNode)j, true);
                logHR += (jOldBranchProb - jNewBranchProb);
            }

            return logHR;
        }

        // randomly selected nodes i and j are not valid candidates for a wide exchange.
        // reject instead of counting (like we do for narrow)
        return Double.NEGATIVE_INFINITY;
    }


    // exchange sub-trees whose root are i and j
    protected void exchangeNodes(Node i, Node j, Node p, Node jP, boolean exchangeEvents) {

        // precondition p -> i & jP -> j
        replace(p, i, j);
        replace(jP, j, i);
        // postcondition p -> j & p -> i

        if (exchangeEvents) {
            double cutHeight = Math.max(i.getHeight(), j.getHeight());
            List<EventNode> iEventsBelow = new ArrayList<>();
            List<EventNode> iEventsAbove = new ArrayList<>();
            for (EventNode e : ((AnnotatedNode)i).getEvents()) {
                if (e.getHeight() <= cutHeight) {
                    iEventsBelow.add(e);
                } else {
                    iEventsAbove.add(e);
                }
            }
            List<EventNode> jEventsBelow = new ArrayList<>();
            List<EventNode> jEventsAbove = new ArrayList<>();
            for (EventNode e : ((AnnotatedNode)j).getEvents()) {
                if (e.getHeight() <= cutHeight) {
                    jEventsBelow.add(e);
                } else {
                    jEventsAbove.add(e);
                }
            }

            ((AnnotatedNode)i).setEvents(iEventsBelow);
            ((AnnotatedNode)i).addEvents(jEventsAbove, true);
            ((AnnotatedNode)j).setEvents(jEventsBelow);
            ((AnnotatedNode)j).addEvents(iEventsAbove, true);
        }
    }

}
