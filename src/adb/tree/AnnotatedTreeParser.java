package adb.tree;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;

import java.util.List;


@Description("Class to initialize an AnnotatedTree from newick tree with type metadata")
public class AnnotatedTreeParser extends AnnotatedTree implements StateNodeInitialiser {

    public Input<String> newickInput =
            new Input<>("newick", "tree in newick format", Input.Validate.REQUIRED);

    public Input<Boolean> adjustTipHeightsInput =
            new Input<>("adjustTipHeights", "adjust tip heights in tree? (default true)", true);

    public Input<Double> scaleInput =
            new Input<>("scale", "factor used to multiply internal node heights during parsing", 1.0);


    public AnnotatedTreeParser() { }


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // read tree
        Tree tree = new TreeParser();
        tree.initByName(
                "newick", newickInput.get(),
                "adjustTipHeights", adjustTipHeightsInput.get(),
                "scale", scaleInput.get(),
                "IsLabelledNewick", true);

        boolean containsEvents = false;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            if (node.getChildCount() == 1) {
                containsEvents = true;
                break;
            }
        }

        if (containsEvents) {
            // tree contains also single-descendant (hidden) nodes
            convertEventTree(tree);
        } else {
            // tree contains only branching nodes and tips
            convertBranchingTree(tree);
        }
    }


    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodeList) {
        stateNodeList.add(this);
    }

}
