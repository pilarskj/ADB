package adbp;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

public class BranchList {
    private List<Branch> branches;

    public BranchList() {
        this.branches = new ArrayList<>();
    }

    public List<Branch> getBranches() {
        return branches;
    }

    // BFS traversal method to add branches to the list
    public void traverseTree(TreeInterface tree, double origin) {
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
            int endNode;
            double endTime;
            if (parent == null) {
                endNode = startNode + 1;
                endTime = origin;
            } else {
                endNode = parent.getNr();
                endTime = parent.getHeight();
            }
            String branchType;
            if (node.isLeaf()) {
                branchType = "external";
            } else {
                branchType = "internal";
            }

            branches.add(new Branch(startNode, endNode, startTime, endTime, branchType));  // Add branch to list

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
    public void assignBranchIndices() {
        int ntips = (int) branches.stream()
                .filter(b -> b.branchType == "external")
                .count();

        int intIndex = 0;
        int extIndex = ntips - 1;
        for (Branch b : branches) {
            if (b.branchType == "internal") {
                b.setIndex(intIndex);
                intIndex++;
            } else if (b.branchType == "external") {
                b.setIndex(extIndex);
                extIndex++;
            }
        }
    }

    // Identify child branches
    public void getChildBranches() {
        for (Branch b : branches) {
            if (b.branchType == "internal") {
                List<Integer> children = branches.stream()
                        .filter(x -> x.endNode == b.startNode) // Filter branches where endNode matches startNode
                        .map(x -> x.branchIndex)               // Get the index of matching branches
                        .toList();
                b.setLeftIndex(children.get(0));
                b.setRightIndex(children.get(1));
            }
        }
    }
}


