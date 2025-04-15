package adbp;

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

    public Input<Double> scaleFactorInput =
            new Input<>("scaleFactor", "scaling factor of shape: ranges from 0 to 1 (close to zero is very large jumps, close to 1.0 is very small jumps)", 0.8);

    public Input<Double> branchProportionInput =
            new Input<>("branchProportion", "proportion of internal branches to be scaled", 0.5);

    // further parameters needed (take current samples)
    public Input<RealParameter> lifetimeInput =
            new Input<>("lifetime", "expected lifetime parameter of the ADB process", Input.Validate.REQUIRED);

    public Input<RealParameter> originInput =
            new Input<>("origin", "time of origin of the process", Input.Validate.REQUIRED);


    Tree tree;
    IntegerParameter shape;
    double branchProportion;
    double scaleFactor;
    RealParameter lifetime;
    RealParameter origin;


    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        shape = shapeInput.get();
        branchProportion = branchProportionInput.get();
        scaleFactor = scaleFactorInput.get();
        lifetime = lifetimeInput.get();
        origin = originInput.get();
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
        int b = (int) (shape.getValue() * scale);
        shape.setValue(b);
        hastingsRatio += -Math.log(scale);

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
        Collections.reverse(internalNodes);
        for (Node n : internalNodes) {
            if (Arrays.stream(randomNodeNumbers).anyMatch(x -> x == n.getNr())) {
                double newBranchLength = gammaDist.sample();
                double newHeight;
                if (n.isRoot()) {
                    newHeight = origin.getValue() - newBranchLength;
                } else {
                    newHeight = n.getParent().getHeight() - newBranchLength;
                }
                // check whether new branch lengths fit in the process
                if (newHeight < 0) {
                    return Double.NEGATIVE_INFINITY;
                }
                n.setHeight(newHeight);
                // based on manipulation, add a value to hastingsRatio
            }
        }

        return hastingsRatio; // Hastings ratio, or -Infinity if rejected
    }
}
