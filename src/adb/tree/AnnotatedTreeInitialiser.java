package adb.tree;

import adb.distribution.Parameterization;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;

import java.util.List;
import java.util.stream.IntStream;

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

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        Tree tree = treeInput.get();
        double scale = scaleInput.get();
        Parameterization parameterization = parameterizationInput.get();

        int nTypes = parameterization.getNTypes();
        Double originTime = parameterization.getOriginTime();
        double lifetime = parameterization.getLifetime(0); // TODO: extend for multi-type
        double shape = parameterization.getShape(0);

        // scale tree
        if (scale != 1.0) {
            for (Node node : tree.getNodesAsArray()) {
                node.setHeight(node.getHeight() * scale);
            }
        }

        // convert to AnnotatedTree
        convertBranchingTree(tree);

        // parsimonous types
        if (nTypes > 1) {
            for (Node node : this.getInternalNodes()) {
                AnnotatedNode left = (AnnotatedNode)node.getLeft();
                AnnotatedNode right = (AnnotatedNode)node.getRight();
                int i = left.getType();
                int j = right.getType();
                int type = 0;
                if (i == j) { // symmetric division
                    for (int k = 1; k < nTypes; i++) { // find most likely parent
                        if (parameterization.getSymTransition(k, i) > parameterization.getSymTransition(type, j)) {
                            type = k;
                        }
                    }
                } else if (i != j) {
                    if (parameterization.getSymTransition(i, j) == 0.0 && parameterization.getSymTransition(j, i) == 0.0) { // impossible transition: add hidden node
                        // TODO
                    } else if (parameterization.getSymTransition(i, j) >= parameterization.getSymTransition(j, i)) {
                        type = i;
                    } else {
                        type = j;
                    }
                }
                AnnotatedNode aNode = (AnnotatedNode)node;
                aNode.setType(type);
            }
        }

        // draw events along branches
        IntStream.range(0, this.getNodeCount())
                .parallel()
                .forEach(i -> {
                    AnnotatedNode node = (AnnotatedNode)this.getNode(i);
                    if (node.isRoot()) {
                        if (originTime != null) {
                            drawEvents(node, node.getHeight(), originTime, lifetime, shape);
                        }
                    } else {
                        drawEvents(node, node.getHeight(), node.getParent().getHeight(), lifetime, shape);
                    }
                });
    }

    // Draw hidden events according to lifetime distribution
    private void drawEvents(AnnotatedNode node, double start, double end, double lifetime, double shape) {
        int type = node.getType();
        node.clearEvents();
        double length = end - start;
        int eventCount = (int)(length / lifetime);
        if (eventCount > 0) {
            double segmentLength = length / (eventCount + 1);
            double time = end - segmentLength;
            while (time > start) {
                node.addEvent(new EventNode(type, time), 1);
                time -= segmentLength;
            }
            node.distributeEvents(shape, end);
        }
    }


    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
}
