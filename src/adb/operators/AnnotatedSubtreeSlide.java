package adb.operators;

import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedNode.EventNode;
import adb.tree.AnnotatedTree;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

// modified version of beast.base.evolution.operator.SubtreeSlide
public class AnnotatedSubtreeSlide extends TreeOperator {

    public Input<Integer> sizeInput = new Input<>("size", "the size of the sliding window up and down, default 1", 1);

    protected int size;

    @Override
    public void initAndValidate() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);
        if (!(tree instanceof AnnotatedTree)) {
            throw new IllegalArgumentException("Attempted to initialise annotated tree operator with regular tree.");
        }

        size = sizeInput.get();
        //limit = limitInput.get();
    }

    @Override
    public double proposal() {

        Tree tree = (Tree) InputUtil.get(treeInput, this);

        double logq;

        Node i;
        final boolean markClades = markCladesInput.get();
        // 1. choose a random node avoiding root
        final int nodeCount = tree.getNodeCount();
        if (nodeCount == 1) {
            // test for degenerate case
            return Double.NEGATIVE_INFINITY;
        }
        do {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot());

        final Node p = i.getParent();
        final Node CiP = getOtherChild(p, i);
        final Node PiP = p.getParent();

        // 2. choose magnitude and direction of the move
        int delta;
        do {
            delta = Randomizer.nextInt(2 * size + 1) - size;
        } while (delta == 0);

        double oldHeight = p.getHeight();
        double newHeight;

        // 3. perform move
        // 3.1 if the move is up
        if (delta > 0) {

            if (delta == ((AnnotatedNode) p).getEventCount()) {
                // no shift possible
                return Double.NEGATIVE_INFINITY;

            // 3.1.1 if topology changes
            } else if (delta > ((AnnotatedNode)p).getEventCount()) {

                if (PiP == null) {
                    // no shift possible
                    return Double.NEGATIVE_INFINITY;
                }

                // find new parent
                Node newParent = PiP;
                Node newChild = p;
                // store events along the path
                List<EventNode> eventPath = new ArrayList<>();
                int x = 0;
                for (int k = 0; k < delta; k++) {
                    eventPath.add(((AnnotatedNode)newChild).getEvent(x));
                    x++;
                    if (x >= ((AnnotatedNode)newChild).getEventCount()) {
                        x = 0; // restore counter
                        newChild = newParent;
                        if (newChild == null) {
                            // no shift possible (no events remain)
                            return Double.NEGATIVE_INFINITY;
                        }
                        if (markClades) newParent.makeDirty(Tree.IS_FILTHY);
                        newParent = newParent.getParent();
                    }
                }

                if (x == 0) {
                    // no shift possible (target event is branching node)
                    return Double.NEGATIVE_INFINITY;
                }

                // get target event
                newHeight = ((AnnotatedNode)newChild).getEvent(x).getHeight();
                p.setHeight(newHeight);

                boolean isRoot = newChild.isRoot();

                // update topology and event lists
                replace(p, CiP, newChild);
                replace(PiP, p, CiP);
                ((AnnotatedNode)i).addEvents(eventPath, true);
                ((AnnotatedNode)CiP).addEvents(((AnnotatedNode)p).getEvents(), true);
                List<EventNode> newParentEvents = ((AnnotatedNode)newChild).getEvents().subList(x, ((AnnotatedNode) newChild).getEventCount());
                ((AnnotatedNode)p).setEvents(newParentEvents);
                newParentEvents.clear();

                // creating a new root
                if (isRoot) {
                    p.setParent(null);
                    tree.setRoot(p);
                }

                // no new root
                else {
                    replace(newParent, newChild, p);
                }

                // count the hypothetical sources of this destination
                int possibleSources = intersectingEdges(newChild, oldHeight, null);
                logq = -Math.log(possibleSources);


            // 3.1.2 if topology does not change
            } else {
                newHeight = ((AnnotatedNode)p).getEvent(delta).getHeight();
                p.setHeight(newHeight);
                // copy events from parent branch to both children
                List<EventNode> pEvents = ((AnnotatedNode)p).getEvents().subList(0, delta);
                ((AnnotatedNode)i).addEvents(pEvents, true);
                ((AnnotatedNode)CiP).addEvents(pEvents, true);
                pEvents.clear();
                logq = 0.0;
            }


        // 3.2 if the move is down
        } else {

            int nEvents = ((AnnotatedNode)i).getEventCount();
            if (Math.abs(delta) >= nEvents) {
                // invalid move (slide below i)
                return Double.NEGATIVE_INFINITY;
            }

            if (Math.abs(delta) == ((AnnotatedNode)CiP).getEventCount()) {
                // no shift possible
                return Double.NEGATIVE_INFINITY;

            // 3.2.1 if topology changes
            } else if (Math.abs(delta) > ((AnnotatedNode)CiP).getEventCount()) {

                // get events along the path
                List<EventNode> eventPath = new ArrayList<>(
                        ((AnnotatedNode)i).getEvents().subList(nEvents + delta, nEvents)
                );

                newHeight = eventPath.get(0).getHeight();

                List<Node> newChildren = new ArrayList<>();
                int possibleDestinations = intersectingEdges(CiP, newHeight, newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                int childIndex = Randomizer.nextInt(newChildren.size());
                Node newChild = newChildren.get(childIndex);
                Node newParent = newChild.getParent();

                // if p was root
                if (p.isRoot()) {
                    // new root is CiP
                    replace(p, CiP, newChild);
                    replace(newParent, newChild, p);

                    CiP.setParent(null);
                    tree.setRoot(CiP);

                } else {
                    replace(p, CiP, newChild);
                    replace(PiP, p, CiP);
                    replace(newParent, newChild, p);
                }

                // event lists
                ((AnnotatedNode)CiP).addEvents(((AnnotatedNode)p).getEvents(), true);
                ((AnnotatedNode)i).getEvents().removeIf(e -> e.getHeight() >= newHeight); // only keep events below cut
                ((AnnotatedNode)newChild).getEvents().removeIf(e -> e.getHeight() >= newHeight);
                eventPath.removeIf(e -> e.getHeight() >= newParent.getHeight()); // only keep events below newParent
                ((AnnotatedNode)p).setEvents(eventPath);

                p.setHeight(newHeight);
                if (markClades) {
                    // make dirty the path from the (down) moved node back up to former parent.
                    Node n = p;
                    while(n != CiP) {
                        n.makeDirty(Tree.IS_FILTHY);
                        n = n.getParent();
                    }
                }

                logq = Math.log(possibleDestinations);


            // 3.2.2 if topology does not change
            } else {
                int ix = ((AnnotatedNode)i).getEventCount() + delta; // delta is negative!
                newHeight =  ((AnnotatedNode)i).getEvents().get(ix).getHeight();
                p.setHeight(newHeight);
                // copy events from slide to parent branch, remove events from children
                List<EventNode> events = ((AnnotatedNode)i).getEvents().subList(ix, ((AnnotatedNode)i).getEventCount());
                ((AnnotatedNode)p).addEvents(events,0);
                events.clear();
                ((AnnotatedNode)CiP).getEvents().removeIf(e -> e.getHeight() >= newHeight);
                logq = 0.0;
            }
        }

        System.out.println(((AnnotatedTree)tree).convertAnnotatedTree(false));
        return logq;
    }


    // copied from regular proposal
    private int intersectingEdges(Node node, double height, List<Node> directChildren) {
        final Node parent = node.getParent();

        if (parent == null) {
            // can happen with non-standard non-mutable trees
            return 0;
        }

        if (parent.getHeight() < height) return 0;

        if (node.getHeight() < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }

        if (node.isLeaf()) {
            return 0;
        } else {
            final int count = intersectingEdges(node.getLeft(), height, directChildren) +
                    intersectingEdges(node.getRight(), height, directChildren);
            return count;
        }
    }

}
