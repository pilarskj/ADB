package adb.operators;

import adb.tree.AnnotatedNode;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;


public class EventSampler extends AnnotatedTreeOperator {

    public Input<Boolean> drawEventCountInput = new Input<>("drawEventCount",
            "draw or fix the number of events along internal branches (default true)", true);

    boolean drawEventCount;


    @Override
    public void initAndValidate() {
        super.initAndValidate();
        drawEventCount = drawEventCountInput.get();
    }


    @Override
    public double proposal() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);
        double oldBranchProb, newBranchProb;
        double logHR = 0;

        // TODO: extend to multi-type version
        for (Node node : tree.getNodesAsArray()) {
            oldBranchProb = getBranchProbability((AnnotatedNode)node, drawEventCount);
            newBranchProb = resampleEvents((AnnotatedNode)node, drawEventCount);
            logHR += (oldBranchProb - newBranchProb);
        }

        // for testing
        //Tree flatTree = ((AnnotatedTree)tree).convertAnnotatedTree(false);
        //System.out.println(flatTree.getRoot().toNewick());
        return logHR;
    }

}
