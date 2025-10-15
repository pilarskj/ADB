package adb;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

import java.util.*;

/*
Class for storing branches of a phylogenetic tree in a list.
Contains functions for traversing a tree starting from the root and recording relevant information for each branch.
Keeps track of start and end time, left and right child, and type of each branch.
 */
public class BranchList {
    private List<Branch> branches;

    public BranchList(TreeInterface tree, double originTime, int originType) {
        branches = new ArrayList<>();
        traverseTree(tree, originTime, originType);
        assignBranchIndices();
        branches.sort(Comparator.comparingInt(b -> b.branchIndex)); // sort by index
        findChildBranches();
    }

    // BFS traversal method to add branches to the list
    private void traverseTree(TreeInterface tree, double originTime, int originType) {
        Node root = tree.getRoot();
        if (root == null) {
            return;
        }

        Queue<Node> queue = new LinkedList<>();
        queue.add(root); // start traversal from the root

        while (!queue.isEmpty()) {
            Node node = queue.poll(); // dequeue current node
            Node parent = node.getParent();

            // get start and end time
            int startNode = node.getNr();
            double startTime = node.getHeight();
            int endNode;
            double endTime;
            if (parent == null) { // stem branch
                endNode = startNode + 1;
                endTime = originTime;
            } else {
                endNode = parent.getNr();
                endTime = parent.getHeight();
            }

            // get start and end type
            int startType;
            int endType;
            Object nodeType = node.getMetaData("type"); // TODO: allow for custom type label
            if (nodeType != null) { startType = ((Double) nodeType).intValue(); } else { startType = 0; } // or -1 // TODO: adapt default types
            if (parent == null) {
                endType = originType;
            } else {
                Object parentType = parent.getMetaData("type");
                if (parentType != null) { endType = ((Double) parentType).intValue(); } else { endType = 0; } // or -1
            }

            // get mode
            String branchMode;
            if (node.isLeaf()) {
                branchMode = "external";
            } else {
                branchMode = "internal";
            }

            branches.add(new Branch(startNode, endNode, startTime, endTime, startType, endType, branchMode)); // add branch to list

            if (node.getLeft() != null) {
                queue.add(node.getLeft()); // enqueue left child
            }
            if (node.getRight() != null) {
                queue.add(node.getRight()); // enqueue right child
            }
        }
    }

    // Assign branch indices
    private void assignBranchIndices() {
        int ntips = countExternalBranches();

        int intIndex = 0;
        int extIndex = ntips - 1;
        for (Branch b : branches) {
            if (b.branchMode.equals("internal")) {
                b.branchIndex = intIndex;
                intIndex++;
            } else if (b.branchMode.equals("external")) {
                b.branchIndex = extIndex;
                extIndex++;
            }
        }
    }

    // Identify child branches
    private void findChildBranches() {
        for (Branch b : branches) {
            if (b.branchMode.equals("internal")) {
                List<Integer> children = branches.stream()
                        .filter(x -> x.endNode == b.startNode) // filter branches where endNode matches startNode
                        .map(x -> x.branchIndex) // get the index of matching branches
                        .toList();
                b.leftIndex = children.get(0);
                b.rightIndex = children.get(1);
            }
        }
    }

    // Get branches as List for looping
    public List<Branch> listBranches() {
        return branches;
    }

    // Get branch with a specific branchIndex
    public Branch getBranchByIndex(int i) {
        for (Branch branch : branches) {
            if (branch.branchIndex == i) {
                return branch;
            }
        }
        throw new IllegalArgumentException("Branch not found for index: " + i);
    }

    // Find branch corresponding to a specific node
    public Branch findBranchForNode(int nodeNr) {
        for (Branch branch : branches) {
            if (branch.startNode == nodeNr) {
                return branch;
            }
        }
        throw new IllegalArgumentException("Branch not found for node: " + nodeNr);
    }

    // Get the number of external branches
    public int countExternalBranches() {
        int i = 0;
        for (Branch branch : branches) {
            if (branch.branchMode.equals("external")) {
                i += 1;
            }
        }
        return i;
    }

    // Get the number of internal branches
    public int countInternalBranches() {
        int i = 0;
        for (Branch branch : branches) {
            if (branch.branchMode.equals("internal")) {
                i += 1;
            }
        }
        return i;
    }

}

