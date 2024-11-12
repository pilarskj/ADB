package adbp;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import org.apache.commons.math3.distribution.GammaDistribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Simulator {
    private Random random;
    private int eventCounter;

    public Simulator() {
        random = new Random();
        random.setSeed(12345); // fixed seed for now - does not seem to work?
        eventCounter = 0;
    }

    // Simulate a phylogeny from ADBP
    public Tree simulateTree(double origin, int originType, double[] a, double[] b, double[] d, double[][] Xsi_as, double[][] Xsi_s, double rho) {
        Tree tree = simulateCompleteTree(origin, originType, a, b, d, Xsi_as, Xsi_s);
        tree = pruneTree(tree, rho);
        return tree;
    }


    // Simulate the full Age-Dependent Branching Process
    public Tree simulateCompleteTree(double origin, int originType, double[] a, double[] b, double[] d, double[][] Xsi_as, double[][] Xsi_s) {

        double rootTime = origin - sampleLifetime(a[originType], b[originType]);

        // Create the root node for the tree (first event)
        Node root = new Node();
        root.setNr(0);  // node number
        root.setHeight(rootTime);  // birth time (height in BEAST2)
        root.setMetaData("type", originType); // type (trait)
        root.metaDataString = "type=" + originType; // add string (for Newick format)

        // Create a Tree object to store the nodes
        Tree tree = new Tree(root);

        // List to keep track of living particles that need to be processed
        List<Node> events = new ArrayList<>();
        events.add(root);

        // Simulate the ADBP process
        while (!events.isEmpty()) {
            Node event = events.remove(0); // look at an event
            int nodeType = (int) event.getMetaData("type"); // type of current particle

            if (random.nextDouble() < 1 - d[nodeType]) { // birth (current particle gives rise to two children)
                // sample types of children
                int[] childrenTypes = sampleChildrenTypes(nodeType, Xsi_as, Xsi_s);
                int leftType = childrenTypes[0];
                int rightType = childrenTypes[1];

                // sample lifetimes and calculate next event times
                double leftHeight = event.getHeight() - sampleLifetime(a[leftType], b[leftType]);
                double rightHeight = event.getHeight() - sampleLifetime(a[rightType], b[rightType]);

                // add children
                Node leftChild = new Node();
                leftChild.setNr(eventCounter++);
                leftChild.setMetaData("type", leftType);
                leftChild.metaDataString = "type=" + leftType;
                event.addChild(leftChild); // add child to node
                tree.addNode(leftChild); // add child to tree
                if (leftHeight < 0) {
                    leftChild.setHeight(0); // this is a tip
                } else {
                    leftChild.setHeight(leftHeight);
                    events.add(leftChild); // add to event list
                }

                Node rightChild = new Node();
                rightChild.setNr(eventCounter++);
                rightChild.setMetaData("type", rightType);
                rightChild.metaDataString = "type=" + rightType;
                event.addChild(rightChild); // add child to node
                tree.addNode(rightChild); // add child to tree
                if (rightHeight < 0) {
                    rightChild.setHeight(0); // this is a tip
                } else {
                    rightChild.setHeight(rightHeight);
                    events.add(rightChild); // add to event list
                }

            }  // in case of death, do nothing
        }

        // Return the constructed tree
        return tree;
    }

    // Prune unsampled and dead particles (and any ancestors without sampled descendants)
    public Tree pruneTree(Tree tree, double rho) {

        List<Node> nodesToPrune = new ArrayList<>();

        // collect leaves to prune
        List<Node> leaves = tree.getExternalNodes();
        for (Node node : leaves) {
            if (node.getHeight() > 0) { // dead particle
                nodesToPrune.add(node);
            } else if (node.getHeight() == 0) {
                if (random.nextDouble() < 1 - rho) { // unsampled particle
                    nodesToPrune.add(node);
                }
            }
        }

        int n = nodesToPrune.size();
        //System.out.println(n + " leaves to prune: ");
        if (n == leaves.size()) {
            //System.out.println("nothing remains!");
            return null;
        }
        /*
        // print the list
        int[] leavesToPrune = new int[n];
        for (int i = 0; i < n; i++) {
            leavesToPrune[i] = nodesToPrune.get(i).getNr();
        }
        System.out.println(Arrays.toString(leavesToPrune));
        */

        // prune ancestors to reconstruct tree structure
        while (!nodesToPrune.isEmpty()) {
            Node node = nodesToPrune.remove(0); // take a node
            if (node.isRoot()) {
                return null;
            }
            Node parent = node.getParent();
            parent.removeChild(node); // deconnect it from parent
            tree.removeNode(node.getNr()); // remove from tree

            if (parent.getChildCount() == 1) {
                Node remainingChild = parent.getChild(0);
                if (parent.isRoot()) {
                    remainingChild.setParent(null);
                    tree.removeNode(parent.getNr());
                    tree.setRoot(remainingChild); // the remaining child becomes the new root
                } else {
                    Node grandParent = parent.getParent();
                    grandParent.removeChild(parent);
                    tree.removeNode(parent.getNr());
                    grandParent.addChild(remainingChild); // reconnect the remaining child to the grandparent
                }
            } else if (parent.getChildCount() == 0) {
                nodesToPrune.add(parent); // continue recursively
            }
            //System.out.println(tree.getNodeCount());
        }

        return tree;
    }


    // Helper method to sample a lifetime from a Gamma distribution
    private double sampleLifetime(double a, double b) {
        GammaDistribution gammaDist = new GammaDistribution(b, a);
        return gammaDist.sample();
    }


    // Helper method to sample offspring types
    private int[] sampleChildrenTypes(int parentType, double[][] Xsi_as, double[][] Xsi_s) {
        // To Do: add matrix checks
        int ntypes = Xsi_s.length;

        double cumulativeProb = 0.0;

        for (int j = 0; j < ntypes; j++) {
            // symmetric case
            cumulativeProb += Xsi_s[parentType][j];
            if (random.nextDouble() < cumulativeProb) {
                return new int[]{j, j};
            }
            // asymmetric case
            cumulativeProb += 2*Xsi_as[parentType][j];
            if (random.nextDouble() < cumulativeProb) {
                if (random.nextDouble() < 0.5) {
                    return new int[]{parentType, j};
                } else {
                    return new int[]{j, parentType};
                }
            }
        }
        return new int[]{-1, -1}; // in case of error
    }
}
