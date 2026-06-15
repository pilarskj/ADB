package adb.operators;

import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedNode.EventNode;
import adb.tree.AnnotatedTree;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

// modified version of beast.base.evolution.operator.WilsonBalding
public class AnnotatedWilsonBalding extends TreeOperator {


    @Override
    public void initAndValidate() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);
        if (!(tree instanceof AnnotatedTree)) {
            throw new IllegalArgumentException("Attempted to initialise annotated tree operator with regular tree.");
        }
    }

    @Override
    public double proposal() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);

        double oldMinAge, newMinAge, newEventHeight, hastingsRatio;
        int oldEventsCount, newEventsCount, oldMinNr, newMinNr, newEventNr;

        // choose a random node avoiding root
        final int nodeCount = tree.getNodeCount();
        Node i;
        do {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot());
        final Node p = i.getParent();

        // choose another random node to insert i above
        Node j;
        Node jP;
        // make sure that the target branch <k, j> is above the subtree being moved
        do {
            j = tree.getNode(Randomizer.nextInt(nodeCount));
            jP = j.getParent();
        } while ((jP != null && jP.getHeight() <= i.getHeight()) || (i.getNr() == j.getNr()));

        // disallow moves that change the root
        if (j.isRoot() || p.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        }

        assert jP != null;  // j != root tested above
        final int pnr = p.getNr();
        final int jPnr = jP.getNr();
        if (jPnr == pnr || j.getNr() == pnr || jPnr == i.getNr()) {
            return Double.NEGATIVE_INFINITY;
        }

        final Node CiP = getOtherChild(p, i);

        // forward move: count available positions
        newMinAge = Math.max(i.getHeight(), j.getHeight());
        newMinNr = ((AnnotatedNode)j).findEventAbove(newMinAge);
        newEventsCount = ((AnnotatedNode)j).getEventCount() - newMinNr;

        // backward move: count available positions (assuming the forward move has happened)
        oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
        oldMinNr = ((AnnotatedNode)CiP).findEventAbove(oldMinAge);
        oldEventsCount = ((AnnotatedNode)CiP).getEventCount() - oldMinNr + ((AnnotatedNode)p).getEventCount();

        if (newEventsCount == 0 || oldEventsCount == 0) {
            // no move
            return Double.NEGATIVE_INFINITY;
        }

        hastingsRatio = (double) newEventsCount / oldEventsCount;

        // sample
        newEventNr = newMinNr + Randomizer.nextInt(newEventsCount);
        newEventHeight = ((AnnotatedNode)j).getEvents().get(newEventNr).getHeight();

        // disconnect p
        final Node pP = p.getParent();
        replace(pP, p, CiP);
        // re-attach, first child node to p
        replace(p, CiP, j);
        // then parent node of j to p
        replace(jP, j, p);

        // mark paths to common ancestor as changed
        if(markCladesInput.get()) {
            Node iup = pP;
            Node jup = p;
            while (iup != jup) {
                if (iup.getHeight() < jup.getHeight()) {
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

        // update height
        p.setHeight(newEventHeight);

        // segregate events
        // collect events below cut (will remain for node j) and above cut (will be transferred to node p)
        List<EventNode> jEventsBelow = new ArrayList<>();
        List<EventNode> jEventsAbove = new ArrayList<>();
        for (EventNode e : ((AnnotatedNode)j).getEvents()) {
            if (e.getHeight() < newEventHeight) {
                jEventsBelow.add(e);
            } else {
                jEventsAbove.add(e);
            }
        }

        ((AnnotatedNode)CiP).addEvents(((AnnotatedNode)p).getEvents(), true); // join lists for CiP and p
        ((AnnotatedNode)i).getEvents().removeIf(e -> e.getHeight() >= newEventHeight); // only keep events below cut
        ((AnnotatedNode)p).setEvents(jEventsAbove);
        ((AnnotatedNode)j).setEvents(jEventsBelow);

        // for testing
        Tree flatTree = ((AnnotatedTree)tree).convertAnnotatedTree(false);
        System.out.println(flatTree.getRoot().toNewick());
        return Math.log(hastingsRatio);
    }

}
