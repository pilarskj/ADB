package adb.tree;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;

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
        initArrays();
    }


    @Override
    protected AnnotatedNode newNode() {
        return new AnnotatedNode();
    }


    // function to convert a strictly bifurcating tree in AnnotatedTree
    // cf. https://github.com/tgvaughan/MultiTypeTree/blob/master/src/multitypetree/evolution/tree/MultiTypeTreeFromUntypedNewick.java
    protected void convertBranchingTree(Tree tree) {

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

            if (node.isRoot()) {
                aNode.setParent(null);
            } else {
                aNode.setParent(annotatedNodes[node.getParent().getNr()]);
            }

            while (aNode.getChildrenMutable().size() < node.getChildCount()) {
                aNode.getChildrenMutable().add(null);
            }

            for (int c = 0; c < node.getChildCount(); c++) {
                aNode.setChild(c, annotatedNodes[node.getChild(c).getNr()]);
            }
        }

        // construct AnnotatedTree
        AnnotatedNode aRoot = annotatedNodes[annotatedNodes.length - 1];
        assignFromWithoutID(new AnnotatedTree(aRoot));
        initArrays();
    }


    // function to convert a tree with single-child nodes in AnnotatedTree
    // cf. https://github.com/tgvaughan/MultiTypeTree/blob/master/src/multitypetree/evolution/tree/MultiTypeTree.java#L538
    // but using tips-to-root traversal
    protected void convertEventTree(Tree tree) {

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

        // construct AnnotatedTree
        assignFromWithoutID(new AnnotatedTree(aRoot));
        initArrays();
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


    // TODO: override Tree functions -- How will those change in BEAST 2.8? Are all these modifications necessary?
    /**
     * ************************ *
     * Methods ported from Tree *
     * ************************ *
     */

    /**
     * Initialise tree-as-array representation + its stored variant *
     */
    @Override
    public void initArrays() {
        // initialise tree-as-array representation + its stored variant
        m_nodes = new AnnotatedNode[nodeCount];
        listNodes((AnnotatedNode)root, (AnnotatedNode[])m_nodes);
        m_storedNodes = new AnnotatedNode[nodeCount];
        Node copy = root.copy();
        listNodes((AnnotatedNode)copy, (AnnotatedNode[])m_storedNodes);
    }

    /**
     * Convert tree to array representation *
     */
    private void listNodes(AnnotatedNode node, AnnotatedNode[] nodes) {
        nodes[node.getNr()] = node;
        node.setTree(this);
        if (!node.isLeaf()) {
            listNodes((AnnotatedNode)node.getLeft(), nodes);
            if (node.getRight()!=null)
                listNodes((AnnotatedNode)node.getRight(), nodes);
        }
    }

    /**
     * Deep copy, returns a completely new tree *
     */
    @Override
    public AnnotatedTree copy() {
        AnnotatedTree tree = new AnnotatedTree();
        tree.ID = ID;
        tree.index = index;
        tree.root = root.copy();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
        return tree;
    }

    /**
     * Copy of all values from existing annotated tree *
     */
    @Override
    public void assignFrom(StateNode other) {
        AnnotatedTree tree = (AnnotatedTree) other;
        AnnotatedNode[] nodes = new AnnotatedNode[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            nodes[i] = new AnnotatedNode();
        }
        ID = tree.ID;
        root = nodes[tree.root.getNr()];
        root.assignFrom(nodes, tree.root);
        root.setParent(null);
        nodeCount = tree.nodeCount;
        internalNodeCount = tree.internalNodeCount;
        leafNodeCount = tree.leafNodeCount;
        initArrays();
    }

    /**
     * As assignFrom, but only copy tree structure *
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        AnnotatedTree tree = (AnnotatedTree) other;
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[tree.root.getNr()];
        Node[] otherNodes = tree.m_nodes;
        int rootNr = root.getNr();
        assignFromFragileHelper(0, rootNr, otherNodes);
        root.setHeight(otherNodes[rootNr].getHeight());
        root.setParent(null);

        AnnotatedNode aRoot = (AnnotatedNode)root;
        aRoot.events.clear();

        if (otherNodes[rootNr].getLeft() != null) {
            root.setLeft(m_nodes[otherNodes[rootNr].getLeft().getNr()]);
        } else {
            root.setLeft(null);
        }
        if (otherNodes[rootNr].getRight() != null) {
            root.setRight(m_nodes[otherNodes[rootNr].getRight().getNr()]);
        } else {
            root.setRight(null);
        }
        assignFromFragileHelper(rootNr + 1, nodeCount, otherNodes);
    }

    /**
     * Helper to assignFromFragile *
     */
    private void assignFromFragileHelper(int start, int end, Node[] otherNodes) {
        for (int i = start; i < end; i++) {
            AnnotatedNode sink = (AnnotatedNode)m_nodes[i];
            AnnotatedNode src = (AnnotatedNode)otherNodes[i];
            sink.setHeight(src.getHeight());
            sink.setParent(m_nodes[src.getParent().getNr()]);
            sink.events.clear();
            sink.events.addAll(src.events);
            if (src.getLeft() != null) {
                sink.setLeft(m_nodes[src.getLeft().getNr()]);
                if (src.getRight() != null) {
                    sink.setRight(m_nodes[src.getRight().getNr()]);
                } else {
                    sink.setRight(null);
                }
            }
        }
    }

    @Override
    protected void store() {
        int rootNr = root.getNr();
        storedRoot = m_storedNodes[rootNr];

        storeNodes(0, rootNr);

        storedRoot.setHeight(m_nodes[rootNr].getHeight());
        storedRoot.setParent(null);
        if (root.getLeft() != null) {
            storedRoot.setLeft(m_storedNodes[root.getLeft().getNr()]);
        } else {
            storedRoot.setLeft(null);
        }
        if (root.getRight() != null) {
            storedRoot.setRight(m_storedNodes[root.getRight().getNr()]);
        } else {
            storedRoot.setRight(null);
        }

        AnnotatedNode aStoredRoot = (AnnotatedNode)storedRoot;
        aStoredRoot.events.clear();
        aStoredRoot.events.addAll(((AnnotatedNode)m_nodes[rootNr]).events);
        storeNodes(rootNr + 1, nodeCount);
    }

    /**
     * Helper to store *
     */
    private void storeNodes(int start, int end) {
        for (int i = start; i<end; i++) {
            AnnotatedNode sink = (AnnotatedNode)m_storedNodes[i];
            AnnotatedNode src = (AnnotatedNode)m_nodes[i];
            sink.setHeight(src.getHeight());
            sink.setParent(m_storedNodes[src.getParent().getNr()]);
            if (src.getLeft() != null) {
                sink.setLeft(m_storedNodes[src.getLeft().getNr()]);
                if (src.getRight() != null)
                    sink.setRight(m_storedNodes[src.getRight().getNr()]);
                else
                    sink.setRight(null);
            }
            sink.events.clear();
            sink.events.addAll(src.events);
        }
    }

    // TODO: String representation of annotated tree

    // TODO: Methods implementing the Loggable interface: init, log, close PrintStream

    // TODO: Function for reconstructing tree from XML fragment in the form of a DOM node

}
