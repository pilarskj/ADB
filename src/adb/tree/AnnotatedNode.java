package adb.tree;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;


@Description("An annotated node in a multi-type phylogenetic tree")
public class AnnotatedNode extends Node {

    // for sorting and comparing events by their timings
    private static final Comparator<EventNode> comparator = (n1, n2) ->
            Double.compare(n1.getHeight(), n2.getHeight());


    // marks all events on the branch above this node (backwards in time, so heights are ascending)
    // this node is first in the list
    protected List<EventNode> events = new ArrayList<>();


    // TODO: necessary?
    public AnnotatedNode() {
        super();
        /* // initialize with this node being the only event
        this.events = new ArrayList<>();
        EventNode branchingEvent = new EventNode(getType(this), this.height);
        this.events.add(branchingEvent); */
    }


    public List<EventNode> getEvents() {
        return events;
    }

    public void setEvents(List<EventNode> events) {
        this.events = events;
    }

    public int getEventCount() {
        return events.size();
    }

    public EventNode getEvent(int idx) {
        return events.get(idx);
    }

    public int getInitialType() { // past type
        return events.get(events.size() - 1).type;
    }

    public int getType() { // present type
        return events.get(0).type;
    }

    // get type of any Node // TODO: make more generic?
    public static int getType(Node node) {
        int type;
        if (node.getMetaData("type") == null) {
            type = 0;
        } else {
            type = (int) (double) node.getMetaData("type");
        }
        return type;
    }

    public void addEvent(EventNode event) {
        int idx = Collections.binarySearch(events, event, comparator);
        if (idx < 0) { idx = -(idx + 1); } // convert to insertion point
        events.add(idx, event);
    }

    public void removeEvent(int idx) {
        if (idx >= events.size())
            throw new IllegalArgumentException("Index to removeEvent() out of range.");
        events.remove(idx);
    }

    public void sortEvents() {
        events.sort(comparator);
    }

    public EventNode getParentEvent(int idx) {
        int pidx = idx + 1;
        if (pidx < events.size()) {
            return events.get(pidx); // predecessor is a hidden node
        } else {
            AnnotatedNode parent = (AnnotatedNode) this.parent;
            if (parent == null) {
                return null; // this event is first in the whole tree
            } else {
                return parent.events.get(0); // predecessor is a branching node
            }
        }
    }

    public List<EventNode> getChildEvents(int idx) { // get next event(s)
        List<EventNode> children = new ArrayList<>();
        if (idx > 0) {
            // one child (hidden node)
            children.add(0, events.get(idx - 1));
        } else {
            if (this.isLeaf()) {
                // no children
                return null;
            } else {
                // two children
                AnnotatedNode left = (AnnotatedNode) this.getLeft();
                List<EventNode> leftEvents = left.getEvents();
                children.add(0, leftEvents.get(leftEvents.size() - 1));

                AnnotatedNode right = (AnnotatedNode) this.getRight();
                List<EventNode> rightEvents = right.getEvents();
                children.add(1, rightEvents.get(rightEvents.size() - 1));
            }
        }
        return children;
    }

    public double[] getWaitingTimes() { // if operating on trees without stem branch
        double[] times = new double[events.size()];
        for (int i = 0; i < times.length; i++) {
            if (this.getParentEvent(i) == null) { // this should never happen
                times[i] = 0;
            } else {
                times[i] = this.getParentEvent(i).getHeight() - events.get(i).getHeight();
            }
        }
        return times;
    }

    public double[] getWaitingTimes(double origin) { // if operating on trees with stem branch
        double[] times = new double[events.size()];
        for (int i = 0; i < times.length; i++) {
            if (this.getParentEvent(i) == null) { // relevant for the first event in the tree
                times[i] = origin - events.get(i).getHeight();
            } else {
                times[i] = this.getParentEvent(i).getHeight() - events.get(i).getHeight();
            }
        }
        return times;
    }


    // TODO: move? change input types to parameters?
    // draw events (without type transitions) along the upstream branch
    public void drawEvents(double lifetime, double shape) {
        events.clear();
        events.add(new EventNode(getType(this), this.getHeight()));

        if (this.isLeaf()) {
            double t = this.getParent().getHeight();
            // work in progress
        }
    }


    /**
     * ************************ *
     * Methods ported from Node *
     * ************************ *
     */

    /**
     * Assign values from a tree in array representation *
     */
    @Override
    public void assignFrom(Node[] nodes, final Node node) {
        height = node.getHeight();
        labelNr = node.getNr();
        metaDataString = node.metaDataString;
        parent = null;
        ID = node.getID();

        AnnotatedNode aNode = (AnnotatedNode)node;
        events.clear();
        events.addAll(aNode.events);

        if (node.getLeft() != null) {
            setLeft(nodes[node.getLeft().getNr()]);
            getLeft().assignFrom(nodes, node.getLeft());
            getLeft().setParent(this);
            if (node.getRight() != null) {
                setRight(nodes[node.getRight().getNr()]);
                getRight().assignFrom(nodes, node.getRight());
                getRight().setParent(this);
            }
        }
    }

    /**
     * (Deep) copy of node *
     */
    public AnnotatedNode copy() {
        AnnotatedNode node = new AnnotatedNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.parent = null;
        node.ID = ID;
        node.events.addAll(events);
        if (getLeft()!=null) {
            node.setLeft(getLeft().copy());
            node.getLeft().setParent(node);
            if (getRight()!=null) {
                node.setRight(getRight().copy());
                node.getRight().setParent(node);
            }
        }
        return node;
    }

}
