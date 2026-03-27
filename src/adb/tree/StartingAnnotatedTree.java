package adb.tree;

import beast.base.core.Description;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;

import java.util.List;

/*
UPGMA + scaling + regular segments + parsimonous type transitions
 */
@Description("Class to initialize an AnnotatedTree")
public class StartingAnnotatedTree extends AnnotatedTree implements StateNodeInitialiser {

    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
}
