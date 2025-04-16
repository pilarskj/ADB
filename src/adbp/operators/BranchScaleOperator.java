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
import java.util.Collections;
import java.util.List;

/*
IDEA:
- joint tree & shape operator
- should make the branch lengths more equal (and thus nodes more synchronous) for high shape, and more diverse for low shape
- branch lengths could be sampled from the lifetime distribution
- resulting trees must be ultrametric, should not exceed the origin, and the topology should not change
- optimizable parameters: scaling factor of shape, proportion of branches updated per move
*/

@Description("Operator that scales internal branch lengths to be more or less equal.")
public class BranchScaleOperator extends Operator {

    public Input<Tree> treeInput =
            new Input<>("tree", "tree", Input.Validate.REQUIRED);

    public Input<IntegerParameter> shapeInput =
            new Input<>("shape", "integer shape parameter of the ADB process", Input.Validate.REQUIRED);
    // will be modified too! add RealParameter option

    final public Input<Integer> windowSizeInput =
            new Input<>("windowSize", "the size of the window both up and down", 1);

//    public Input<Double> scaleFactorInput =
//            new Input<>("scaleFactor", "scaling factor of shape: ranges from 0 to 1 (close to zero is very large jumps, close to 1.0 is very small jumps)", 0.8);

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
    double scaleFactor;
    int windowSize;
    RealParameter lifetime;
    RealParameter origin;


    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        shape = shapeInput.get();
        branchProportion = branchProportionInput.get();
//        scaleFactor = scaleFactorInput.get();
        lifetime = lifetimeInput.get();
        origin = originInput.get();
        windowSize = windowSizeInput.get();
    }

    // cf. https://github.com/CompEvol/beast2/blob/master/src/beast/base/evolution/operator/ScaleOperator.java
    protected double getScaler() {
        return (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
    }


    @Override
    public double proposal() {

        double hastingsRatio = 0;

        // new shape: use scaling or another move?
        double scale = getScaler();
//        int b = (int) (shape.getValue() * scale);
        int b = shape.getValue() + Randomizer.nextInt(2 * windowSize + 1) - windowSize;
        if (b < 1) {
            // invalid move, can be rejected immediately
            return Double.NEGATIVE_INFINITY;
        }
        shape.setValue(b);

//        hastingsRatio -= Math.log(scale);

        // define branch length distribution
        double a = lifetime.getValue() / b;
        GammaDistribution gammaDist = new GammaDistribution(b, a);

        // randomly select internal nodes, operate on branches connecting them to parents
        int branchCount = (int) (branchProportion * (tree.getInternalNodeCount() - 1));
        List<Node> internalNodes = tree.getInternalNodes();
        int[] internalNodeNumbers = internalNodes.stream().mapToInt(Node::getNr).toArray();
        Randomizer.shuffle(internalNodeNumbers);
        int[] randomNodeNumbers = Arrays.copyOfRange(internalNodeNumbers, 0, branchCount);

        // traverse nodes starting from the root, update selected branches
//        Collections.reverse(internalNodes); // you can instead avoid reversing by using the index, saves time
        for (int i = internalNodes.size(); i > 0; i--) {
            Node n = internalNodes.get(i - 1);
            if (Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getNr())) {
                double oldBranchLength = n.getLength();
                double newBranchLength = gammaDist.sample();
                double newHeight;
                if (n.isRoot()) {
                    newHeight = origin.getValue() - newBranchLength;
                } else {
                    newHeight = n.getParent().getHeight() - newBranchLength;
                }
                // check whether new branch lengths fit in the process
                if (newHeight < 0 || checkForChildrenNotOnList(randomNodeNumbers, n)) {
                    return Double.NEGATIVE_INFINITY;
                }
                n.setHeight(newHeight);
                // based on manipulation, add a value to hastingsRatio
                double toCheck = Math.log(gammaDist.density(oldBranchLength) / gammaDist.density(newBranchLength));
                double toCheck2 = (b-1) * (Math.log(oldBranchLength) - Math.log(newBranchLength)) + (newBranchLength-oldBranchLength)/a;
                hastingsRatio += toCheck;
            }
        }

        return hastingsRatio; // Hastings ratio, or -Infinity if rejected
    }

    private boolean checkForChildrenNotOnList(int[] randomNodeNumbers, Node n){
        if (!Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getLeft().getNr()) && (n.getLeft().getHeight() > n.getHeight()))
            return true;
        if (!Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getRight().getNr()) && (n.getRight().getHeight() > n.getHeight()))
            return true;
        return false;
    }

}
