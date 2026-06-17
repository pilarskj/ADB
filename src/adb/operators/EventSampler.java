package adb.operators;

import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedTree;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;


public class EventSampler extends AnnotatedTreeOperator {

    @Override
    public double proposal() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);
        double oldBranchProb, newBranchProb;
        double logHR = 0;

        // TODO: extend to multi-type version
        for (Node node : tree.getNodesAsArray()) {
        //for (Node node : tree.getExternalNodes()) {
            oldBranchProb = getBranchProbability((AnnotatedNode)node);
            newBranchProb = resampleEvents((AnnotatedNode)node);
            logHR += (newBranchProb - oldBranchProb);
        }

        // for testing
        //Tree flatTree = ((AnnotatedTree)tree).convertAnnotatedTree(false);
        //System.out.println(flatTree.getRoot().toNewick());
        return logHR;
    }

}
