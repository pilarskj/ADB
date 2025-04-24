package adbp.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import org.apache.commons.math3.distribution.GammaDistribution;

import java.util.Arrays;
import java.util.List;

/*
IDEA:
- joint tree & shape operator
- makes the branch lengths more equal (and thus nodes more synchronous) for high shape, and more diverse for low shape
- branch lengths are sampled directly from the lifetime distribution
- the topology of the resulting trees does not change
- optimizable parameters: window size for shape proposal, proportion of branches updated per move
*/

@Description("Operator that proposes a new shape parameter and" +
        "samples internal branch lengths from lifetime distribution.")
public class SampleBranchesOperator extends Operator {

    public Input<Tree> treeInput =
            new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public Input<IntegerParameter> shapeInput =
            new Input<>("shape", "integer shape parameter of the ADB process", Input.Validate.REQUIRED);
    // TODO add RealParameter option

    final public Input<Integer> windowSizeInput =
            new Input<>("windowSize", "the size of the window both up and down", 5);

    // public Input<Double> scaleFactorInput =
    //         new Input<>("scaleFactor", "scaling factor of shape: ranges from 0 to 1 (close to zero is very large jumps, close to 1.0 is very small jumps)", 0.8);

    public Input<Double> branchProportionInput =
            new Input<>("branchProportion", "proportion of internal branches to be scaled", 0.1);

    // further parameters needed (take current samples)
    public Input<RealParameter> lifetimeInput =
            new Input<>("lifetime", "expected lifetime parameter of the ADB process", Input.Validate.REQUIRED);

    public Input<RealParameter> originInput =
            new Input<>("origin", "time of origin of the process", Input.Validate.REQUIRED);


    Tree tree;
    IntegerParameter shape;
    double branchProportion;
    // double scaleFactor;
    int windowSize;
    RealParameter lifetime;
    RealParameter origin;


    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        shape = shapeInput.get();
        branchProportion = branchProportionInput.get();
        // scaleFactor = scaleFactorInput.get();
        lifetime = lifetimeInput.get();
        origin = originInput.get();
        windowSize = windowSizeInput.get();
    }

    /* // cf. https://github.com/CompEvol/beast2/blob/master/src/beast/base/evolution/operator/ScaleOperator.java
    protected double getScaler() {
        return (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
    } */


    @Override
    public double proposal() {

        double hastingsRatio = 0;

        // old params
        double C = lifetime.getValue();
        int oldShape = shape.getValue();
        GammaDistribution oldGammaDist = new GammaDistribution(oldShape, C/oldShape);

        // new shape
        // double scale = getScaler();
        // int newShape = (int) (b * scale);
        int newShape = oldShape + Randomizer.nextInt(2 * windowSize + 1) - windowSize; // use random walk
        if (newShape < shape.getLower() || newShape > shape.getUpper()) {
            // invalid move, can be rejected immediately
            return Double.NEGATIVE_INFINITY;
        }
        if (newShape == oldShape) {
            // this saves calculating the posterior
            return Double.NEGATIVE_INFINITY;
        }
        // hastingsRatio += Math.log(scale); for RandomWalk, hR = 0
        shape.setValue(newShape);
        // define new branch length distribution
        GammaDistribution newGammaDist = new GammaDistribution(newShape, C/newShape);

        // randomly select internal nodes, operate on branches connecting them to parents
        int branchCount = (int) (branchProportion * (tree.getInternalNodeCount() - 1));
        List<Node> internalNodes = tree.getInternalNodes();
        int[] internalNodeNumbers = internalNodes.stream().mapToInt(Node::getNr).toArray();
        Randomizer.shuffle(internalNodeNumbers);
        int[] randomNodeNumbers = Arrays.copyOfRange(internalNodeNumbers, 0, branchCount);

        /* // alternative: randomly select an internal node, modify branch lengths in its subtree (min cherry, max whole tree)
        // cf. https://github.com/CompEvol/beast2/blob/master/src/beast/base/evolution/operator/Uniform.java
        final int nodeCount = tree.getNodeCount();
        Node node;
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isLeaf());
        List<Node> subtree = node.getAllChildNodesAndSelf(); // recursive, so list is already sorted from node to tips */

        // traverse nodes starting from the root, update selected branches
        // for (Node n : subtree) {
            // if (!n.isLeaf()) {
        for (int i = internalNodes.size(); i > 0; i--) {
            Node n = internalNodes.get(i - 1);
            if (Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getNr())) {
                double oldBranchLength;
                double newBranchLength = newGammaDist.sample();
                double newHeight;
                if (n.isRoot()) {
                    oldBranchLength = origin.getValue() - n.getHeight();
                    newHeight = origin.getValue() - newBranchLength;
                } else {
                    oldBranchLength = n.getLength();
                    newHeight = n.getParent().getHeight() - newBranchLength;
                }
                n.setHeight(newHeight);
                // based on manipulation, add a value to hastingsRatio
                double toCheck = Math.log(newGammaDist.density(newBranchLength)) - Math.log(oldGammaDist.density(oldBranchLength)); // TODO or other way around?
                //double toCheck2 = (b-1) * (Math.log(oldBranchLength) - Math.log(newBranchLength)) + (newBranchLength-oldBranchLength)/a;
                hastingsRatio += toCheck;
            }
        }
        // post-hoc, check violations in tree topology (otherwise, newHeights of children might still be higher than parent)
        // alternative: traverse the whole tree and check for negative lengths
        // TODO isn't there a better way?
        for (Node n : tree.getInternalNodes()) {
            // if (n.getLength() < 0) {
            if (n.getHeight() < 0 || n.getHeight() < n.getLeft().getHeight() || n.getHeight() < n.getRight().getHeight()) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        return hastingsRatio;
    }
}
