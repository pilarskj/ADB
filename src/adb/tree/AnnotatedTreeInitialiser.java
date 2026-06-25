package adb.tree;

import adb.tree.AnnotatedNode.EventNode;
import adb.distribution.Parameterization;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

import static adb.tree.AnnotatedNode.distributeEvents;

/*
UPGMA + scaling + regular segments + parsimonous type transitions
 */
@Description("Class to initialize an AnnotatedTree from Tree")
public class AnnotatedTreeInitialiser extends AnnotatedTree implements StateNodeInitialiser {

    public Input<Tree> branchingTreeInput =
            new Input<>("branchingTree", "a standard BEAST2 branching tree", Input.Validate.REQUIRED);

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);

    public Input<Double> scaleInput =
            new Input<>("scale", "factor used to multiply internal node heights", 1.0);
    // TODO: potentially use automatic scaling to adapt branch lengths to lifetimes


    Parameterization parameterization;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        Tree branchingTree = branchingTreeInput.get();
        double scale = scaleInput.get();
        parameterization = parameterizationInput.get();

        // scale tree
        if (scale != 1.0) {
            for (Node node : branchingTree.getNodesAsArray()) {
                node.setHeight(node.getHeight() * scale);
            }
        }

        // convert to AnnotatedTree
        convertBranchingTree(branchingTree);

        // parsimonous ancestral types
        if (parameterization.getNTypes() > 1) {
            for (Node node : getInternalNodes()) {
                addAncestralTypes((AnnotatedNode)node);
            }
        }

        // draw events along branches
        /* IntStream.range(0, getNodeCount())
                .parallel()
                .forEach(i -> { */
        for (int i = 0; i < getNodeCount(); i++) {
            AnnotatedNode node = (AnnotatedNode)getNode(i);
            if (node.isRoot()) {
                if (parameterization.getOriginTime() != null) {
                    drawEvents(node, parameterization.getOriginTime());
                }
            } else {
                drawEvents(node, node.getParent().getHeight());
            }
            //});
        }

        if (!parameterization.isDirectProgenitor(
                parameterization.getOriginType(),
                ((AnnotatedNode)getRoot()).getInitialType())) {
            Log.warning("WARNING: Type at origin is incompatible with the starting tree.");
        }
    }


    // Draw events along the (multi-type) branch upstream of a node
    private void drawEvents(AnnotatedNode node, double origin) {
        int type = node.getType();
        double start = node.getHeight();
        double end = start;
        int nEvents = node.getEventCount();
        int ix = 1;
        int i = 0;
        while (i < node.getEventCount()) {
            if (node.getEvent(i).getType() == type) {
                if (node.getParentEvent(i) != null) {
                    end = node.getParentEvent(i).getHeight();
                } else {
                    end = origin;
                }
                i++;
            } else {
                // draw for current type
                drawEvents(node, ix, start, end, type);
                i = node.getEventCount() - nEvents + 1;
                // reset start, type, and counters
                start = node.getEvent(i).getHeight();
                type = node.getEvent(i).getType();
                nEvents = node.getEventCount();
                ix += i;
            }
        }
        drawEvents(node, ix, start, end, type);
    }


    // Draw hidden events along a (single-type) branch segment according to lifetime distribution
    private void drawEvents(AnnotatedNode node, int ix, double start, double end, int type) {
        double lifetime = parameterization.getLifetime(type);
        double shape = parameterization.getShape(type);
        double length = end - start;
        int eventCount = (int)(length / lifetime);
        if (eventCount > 0) {
            List<EventNode> events = new ArrayList<>();
            double segmentLength = length / (eventCount + 1);
            double time = end - segmentLength;
            while (time > start) {
                events.add(0, new EventNode(type, time));
                time -= segmentLength;
            }

            distributeEvents(events, start, end, shape);
            for (int i = 0; i < eventCount; i++) {
                node.addEvent(events.get(i), ix);
                ix++;
            }
        }
    }


    // Propagate types from tips to root, adding minimal intermediate transitions
    private void addAncestralTypes(AnnotatedNode node) {
        int nTypes = parameterization.getNTypes();

        int type;
        while (true) {
            List<EventNode> children = node.getChildEvents(0);
            int i = children.get(0).getType();
            int j = children.get(1).getType();

            // determine node type
            if (i == j) { // symmetric division
                type = 0;
                for (int k = 0; k < nTypes; k++) {
                    if (parameterization.getSymTransition(k, i) > parameterization.getSymTransition(type, i)) {
                        type = k;
                    }
                }
                break;
            } else { // asymmetric division
                // resolve impossible transitions
                if (parameterization.getAsymTransition(i, j) == 0.0 && parameterization.getAsymTransition(j, i) == 0.0) {
                    // select downstream branch to add intermediate types
                    int child = parameterization.isProgenitor(i, j) ? 1 : 0;

                    // sample progenitor and height
                    ArrayList<Integer> progenitors = parameterization.getDirectProgenitors(children.get(child).getType());
                    int k;
                    do {
                        k = progenitors.get(Randomizer.nextInt(progenitors.size()));
                    } while (k == children.get(child).getType());
                    double height = Randomizer.uniform(children.get(child).getHeight(), node.getHeight());
                    EventNode event = new EventNode(k, height);
                    ((AnnotatedNode)(node.getChild(child))).addEvent(event, true);

                } else {
                    type = (parameterization.getAsymTransition(i, j) >= parameterization.getAsymTransition(j, i)) ? i : j;
                    break;
                }
            }
        }
        node.setType(type);
    }


    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodeList) {
        stateNodeList.add(this);
    }
}
