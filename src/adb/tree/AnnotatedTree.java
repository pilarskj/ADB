package adb.tree;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.*;

import static adb.tree.AnnotatedNode.getType;


@Description("An annotated multi-type phylogenetic tree.")
public class AnnotatedTree extends Tree {

    public AnnotatedTree() { };

    public AnnotatedTree(Node root) {

        if (!(root instanceof AnnotatedNode))
            throw new IllegalArgumentException("Attempted to instantiate "
                    + "annotated tree with regular root node.");

        setRoot(root); // also updates nodeCount
        internalNodeCount = root.getInternalNodeCount();
        leafNodeCount = root.getLeafNodeCount();
        initArrays(); // TODO: check difference between m_nodes and m_storedNodes
    }


    @Override
    protected AnnotatedNode newNode() {
        return new AnnotatedNode();
    }


    // function to convert a strictly bifurcating tree in AnnotatedTree
    // cf. https://github.com/tgvaughan/MultiTypeTree/blob/master/src/multitypetree/evolution/tree/MultiTypeTreeFromUntypedNewick.java
    protected AnnotatedTree convertBranchingTree(Tree tree) { // TODO: void! assignFrom/ copy?

        // create all nodes
        AnnotatedNode[] annotatedNodes = new AnnotatedNode[tree.getNodeCount()];
        for (int i = 0; i < annotatedNodes.length; i++) {
            Node node = tree.getNode(i);
            AnnotatedNode aNode = new AnnotatedNode();
            aNode.setNr(i);
            aNode.setHeight(node.getHeight());
            aNode.setID(node.getID());
            List<EventNode> events = new ArrayList<>();
            EventNode branchingEvent = new EventNode(getType(node), node.getHeight());
            events.add(branchingEvent);
            aNode.setEvents(events);
            annotatedNodes[i] = aNode;
        }

        // preserve branching structure
        for (int i = 0; i < annotatedNodes.length; i++) {
            AnnotatedNode aNode = annotatedNodes[i];
            Node node = tree.getNode(i);

            if (node.isRoot())
                aNode.setParent(null);
            else
                aNode.setParent(annotatedNodes[node.getParent().getNr()]);

            while (aNode.getChildrenMutable().size() < node.getChildCount()) {
                aNode.getChildrenMutable().add(null);
            }

            for (int c = 0; c < node.getChildCount(); c++) {
                aNode.setChild(c, annotatedNodes[node.getChild(c).getNr()]);
            }
        }

        // construct AnnotatedTree
        AnnotatedNode aRoot = annotatedNodes[annotatedNodes.length-1];
        return new AnnotatedTree(aRoot);

        // if in-place: assign tree topology
        // assignFromWithoutID(new AnnotatedTree(aRoot));
        // initArrays();
    }


    // function to convert a tree with single-child nodes in AnnotatedTree
    // cf. https://github.com/tgvaughan/MultiTypeTree/blob/master/src/multitypetree/evolution/tree/MultiTypeTree.java#L538
    // but using tips-to-root traversal
    protected AnnotatedTree convertEventTree(Tree tree) {
        // map to keep track of the AnnotatedNodes (key: original Node Nr, value: new AnnotatedNode)
        Map<Integer, AnnotatedNode> nodeMap = new HashMap<>();

        // get nodes sorted by height (tips first)
        Node[] nodes = tree.getNodesAsArray();

        // initialise counter of nodes
        AnnotatedNode aRoot = null;
        int nextNr = 0;

        for (Node node : nodes) {
            // only create AnnotatedNode for tips or branching points
            if (node.getChildCount() == 1) continue;

            AnnotatedNode aNode = new AnnotatedNode();
            aNode.setHeight(node.getHeight());
            aNode.setID(node.getID());
            aNode.setNr(nextNr);

            List<EventNode> events = new ArrayList<>();
            events.add(new EventNode(getType(node), node.getHeight()));
            // CLIMB UP: look at the lineage above this node and collect all single-child nodes as events for THIS node
            Node parent = node.getParent();
            while (parent != null && parent.getChildCount() == 1) {
                events.add(new EventNode(getType(parent), parent.getHeight()));
                parent = parent.getParent();
            }
            aNode.setEvents(events);

            // BRANCHING: store this node so its parent can find it later
            if (node.getChildCount() == 2) {
                // find the children we already created and add them to this parent
                aNode.addChild(findBranchingChild(node.getLeft(), nodeMap));
                aNode.addChild(findBranchingChild(node.getRight(), nodeMap));
            }

            nodeMap.put(node.getNr(), aNode);
            aRoot = aNode;
            nextNr++;
        }

        return new AnnotatedTree(aRoot);
    }

    /**
     * Helper to skip the single-child event nodes and find the underlying AnnotatedNode child stored in the map.
     */
    private AnnotatedNode findBranchingChild(Node child, Map<Integer, AnnotatedNode> nodeMap) {
        Node current = child;
        while (current.getChildCount() == 1) {
            current = current.getLeft();
        }
        return nodeMap.get(current.getNr());
    }

}
