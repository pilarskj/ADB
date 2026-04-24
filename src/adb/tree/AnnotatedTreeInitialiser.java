package adb.tree;

import adb.tree.AnnotatedNode.EventNode;
import adb.distribution.Parameterization;
import adb.distribution.Parameterization.TypeMap;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import static adb.tree.AnnotatedNode.distributeEvents;

/*
UPGMA + scaling + regular segments + parsimonous type transitions
 */
@Description("Class to initialize an AnnotatedTree from Tree")
public class AnnotatedTreeInitialiser extends AnnotatedTree implements StateNodeInitialiser {

    public Input<Tree> treeInput =
            new Input<>("tree", "a standard BEAST2 branching tree", Input.Validate.REQUIRED);

    public Input<Double> scaleInput =
            new Input<>("scale", "factor used to multiply internal node heights", 1.0);

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);

    int nTypes;
    Parameterization parameterization;
    TypeMap[] typeMap;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        Tree tree = treeInput.get();
        double scale = scaleInput.get();
        parameterization = parameterizationInput.get();

        int nTypes = parameterization.getNTypes();
        Double originTime = parameterization.getOriginTime();
        typeMap = parameterization.getTypeMap();

        // scale tree
        if (scale != 1.0) {
            for (Node node : tree.getNodesAsArray()) {
                node.setHeight(node.getHeight() * scale);
            }
        }

        // convert to AnnotatedTree
        convertBranchingTree(tree);

        // TODO: solve bug: negative branch lengths in multi-type case
        // parsimonous types
        if (nTypes > 1) {
            for (Node node : getInternalNodes()) {
                addAncestralTypes((AnnotatedNode)node);
            }
        }

        // draw events along branches
        /* IntStream.range(0, getNodeCount())
                .parallel()
                .forEach(i -> { */
        for (int i = 0; i < getNodeCount(); i++) {
            AnnotatedNode node = (AnnotatedNode) getNode(i);
            if (node.isRoot()) {
                if (originTime != null) {
                    drawEvents(node, originTime);
                }
            } else {
                drawEvents(node, node.getParent().getHeight());
            }
            //});
        }
    }


    private void drawEvents(AnnotatedNode node, double origin) {
        int type = node.getType();
        double start = node.getHeight();
        double end = start;
        for (int i = 0; i < node.getEventCount(); i++) {
            if (node.getEvent(i).getType() == type) {
                if (node.getParentEvent(i) != null) {
                    end = node.getParentEvent(i).getHeight();
                } else {
                    end = origin;
                }
            } else {
                // draw for current type
                drawEvents(node, start, end, type);
                // reset start and type
                start = node.getEvent(i).getHeight();
                type = node.getEvent(i).getType();
            }
        }
        drawEvents(node, start, end, type);
    }


    // Draw hidden events according to lifetime distribution
    private void drawEvents(AnnotatedNode node, double start, double end, int type) {
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
                node.addEvent(events.get(i), true); // TODO: should append at provided index ix?
            }
        }
    }

    private void addAncestralTypes(AnnotatedNode node) {
        List<EventNode> children = node.getChildEvents(0);
        int i = children.get(0).getType();
        int j = children.get(1).getType();
        int type = 0;
        if (i == j) { // symmetric division
            for (int k = 1; k < nTypes; i++) { // find most likely parent
                if (parameterization.getSymTransition(k, i) > parameterization.getSymTransition(type, j)) {
                    type = k;
                }
            }
        } else if (i != j) {
            if (parameterization.getSymTransition(i, j) == 0.0 && parameterization.getSymTransition(j, i) == 0.0) { // impossible transition: add hidden node
                // randomly select child, find progenitor type, and add event node
                int child = Randomizer.nextInt(2); // 0 or 1, i.e. left or right
                int k = 0;
                while (!typeMap[k].hasDescendant(j)) { k++; } // TODO: should consider i or j --> StackOverFlow
                EventNode event = new EventNode(k, Randomizer.uniform(children.get(child).getHeight(), node.getHeight()));
                ((AnnotatedNode)(node.getChild(child))).addEvent(event, true);
                // iterate until valid transitions found
                addAncestralTypes(node);
            } else if (parameterization.getSymTransition(i, j) >= parameterization.getSymTransition(j, i)) {
                type = i;
            } else {
                type = j;
            }
        }
        node.setType(type);
    }


    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
}
