package adb.operators;

import adb.tree.AnnotatedNode;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;


public class EventSampler extends AnnotatedTreeOperator {

    public Input<Boolean> drawEventCountInput = new Input<>("drawEventCount",
            "draw or fix the number of events along internal branches (default true)", true);

    public Input<Double> proportionBranchesInput = new Input<>("proportionBranches",
            "proportion of branches to sample events on (default 0.1)", 0.1);

    boolean drawEventCount;
    double proportionBranches;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        drawEventCount = drawEventCountInput.get();
        proportionBranches = proportionBranchesInput.get();
    }


    @Override
    public double proposal() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);
        double oldBranchProb, newBranchProb;
        double logHR = 0;

        double u;
        // TODO: extend to multi-type version
        for (Node node : tree.getNodesAsArray()) {
            u = Randomizer.nextDouble();
            if (u < proportionBranches) {
                oldBranchProb = getBranchProbability((AnnotatedNode)node, drawEventCount);
                newBranchProb = resampleEvents((AnnotatedNode)node, drawEventCount);
                logHR += (oldBranchProb - newBranchProb);
            }
        }

        // for testing
        //Tree flatTree = ((AnnotatedTree)tree).convertAnnotatedTree(false);
        //System.out.println(flatTree.getRoot().toNewick());
        return logHR;
    }

}
