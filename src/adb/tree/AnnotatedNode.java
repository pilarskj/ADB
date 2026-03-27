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

    public int getType() {
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
            throw new IllegalArgumentException("Index to removeChange() out of range.");
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

}
