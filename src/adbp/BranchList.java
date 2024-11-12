package adbp;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;

public class BranchList {
    private List<Branch> branches;

    public BranchList(TreeInterface tree, double origin, boolean typed) {
        branches = new ArrayList<>();
        traverseTree(tree, origin, typed);
        assignBranchIndices();
        findChildBranches();
    }

    // BFS traversal method to add branches to the list
    private void traverseTree(TreeInterface tree, double origin, boolean typed) {
        Node root = tree.getRoot();
        if (root == null) {
            return;
        }

        Queue<Node> queue = new LinkedList<>();
        queue.add(root);  // Start the traversal from the root

        while (!queue.isEmpty()) {
            Node node = queue.poll();  // Dequeue the current node
            Node parent = node.getParent();

            int startNode = node.getNr();
            double startTime = node.getHeight();

            // extract types of a typed tree
            int nodeType;
            if (typed) {
                Double nodeTypeDouble = (Double) node.getMetaData("type");
                nodeType = nodeTypeDouble.intValue();
            } else {
                nodeType = -1;
            }

            int endNode;
            double endTime;
            if (parent == null) {
                endNode = startNode + 1;
                endTime = origin;
            } else {
                endNode = parent.getNr();
                endTime = parent.getHeight();
            }
            String branchMode;
            if (node.isLeaf()) {
                branchMode = "external";
            } else {
                branchMode = "internal";
            }

            branches.add(new Branch(startNode, endNode, startTime, endTime, branchMode, nodeType));  // Add branch to list

            // Process left child
            if (node.getLeft() != null) {
                queue.add(node.getLeft());  // Enqueue the left child
            }

            // Process right child
            if (node.getRight() != null) {
                queue.add(node.getRight());  // Enqueue the right child
            }
        }
    }

    // Assign branch indices
    private void assignBranchIndices() {
        int ntips = (int) branches.stream()
                .filter(b -> b.branchMode.equals("external"))
                .count();

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
                        .filter(x -> x.endNode == b.startNode) // Filter branches where endNode matches startNode
                        .map(x -> x.branchIndex)               // Get the index of matching branches
                        .toList();
                b.leftIndex = children.get(0);
                b.rightIndex = children.get(1);
            }
        }
    }

    // Get branches as List or looping
    public List<Branch> listBranches() {
        return branches;
    }

    // Get branch with a specific branchIndex
    public  Branch getBranchByIndex(int i) {
        for (Branch branch : branches) {
            if (branch.branchIndex == i) {
                return branch;
            }
        }
        // If no branch is found with the given branchIndex, throw an exception
        throw new IllegalArgumentException("Branch not found for index: " + i);
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


    // Assign random types (for testing)
    public void assignRandomTypes(int ntypes) {
        Random random = new Random();
        random.setSeed(12345);
        for (Branch branch : branches) {
            int randomType = random.nextInt(ntypes);
            branch.nodeType = randomType;
        }
    }

}


