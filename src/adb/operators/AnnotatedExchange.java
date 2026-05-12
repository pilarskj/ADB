package adb.operators;

import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedNode.EventNode;
import adb.tree.AnnotatedTree;
import beast.base.evolution.operator.Exchange;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;

import java.util.ArrayList;
import java.util.List;

// modified version of beast.base.evolution.operator.Exchange
public class AnnotatedExchange extends Exchange {

    @Override
    public void initAndValidate() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);
        if (!(tree instanceof AnnotatedTree)) {
            throw new IllegalArgumentException("Attempted to initialise annotated tree operator with regular tree.");
        }
    }

    @Override
    /* exchange sub-trees whose root are i and j */
    protected void exchangeNodes(Node i, Node j,
                                 Node p, Node jP) {
        // precondition p -> i & jP -> j
        replace(p, i, j);
        replace(jP, j, i);
        // postcondition p -> j & p -> i

        // exchange events
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
