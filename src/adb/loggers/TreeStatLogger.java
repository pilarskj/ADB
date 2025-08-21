package adbp.loggers;


import beast.base.core.*;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;


@Description("Logger to report statistics of a tree")
public class TreeStatLogger extends CalculationNode implements Loggable, Function {
    final public Input<Tree> treeInput = new Input<>("tree", "tree to report height for.", Validate.REQUIRED);
    final public Input<Boolean> logHeightInput = new Input<>("logHeight", "If true, tree height will be logged.", true);
    final public Input<Boolean> logLengthInput = new Input<>("logLength", "If true, tree length will be logged.", true);
    final public Input<Boolean> logBranchLengthDistrInput = new Input<>("logBranchLengthDistr",
            "If true, mean and sd for branch lengths will be logged. It excludes the branches leading to tips", true);

    boolean logHeight, logLength, logBranchLengthDistr;

    @Override
    public void initAndValidate() {

        // This confusing line is because the default situation is to log
        // the tree height: we want the height to be logged only if neither of
        // these defaults is altered.
        logHeight = logHeightInput.get();
        logLength = logLengthInput.get();
        logBranchLengthDistr = logBranchLengthDistrInput.get();

    	if (!logHeight && !logLength && !logBranchLengthDistr) {
    		Log.warning.println("TreeStatLogger " + getID() + "logs nothing. Set logHeight=true or logLength=true or logVariance=true to log at least something");
    	}
    }

    @Override
    public void init(PrintStream out) {
        final Tree tree = treeInput.get();
        if (logHeight) {
            out.print(tree.getID() + ".height\t");
        }
        if (logLength) {
            out.print(tree.getID() + ".treeLength\t");
        }
        if (logBranchLengthDistr) {
            out.print(tree.getID() + ".branchLengthMean\t");
            out.print(tree.getID() + ".branchLengthSD\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        final Tree tree = treeInput.get();
        if (logHeight) {
        	out.print(tree.getRoot().getHeight() + "\t");
        }
        if (logLength || logBranchLengthDistr){
            double length = getLength(tree);
            if (logLength) {
                out.print(length + "\t");
            }
            if (logBranchLengthDistr) {
                double internalLength = length - getLeafEdgeLength(tree);
                out.print(internalLength / (tree.getInternalNodeCount() - 1) + "\t");
                double variance = getBranchLengthSD(tree, internalLength);
                out.print(variance + "\t");
            }
        }
    }

    private double getLength(Tree tree) {
    	double length = 0;
    	for (Node node : tree.getNodesAsArray()) {
    		if (!node.isRoot()) {
    			length += node.getLength();
    		}
    	}
		return length;
	}

    private double getLeafEdgeLength(Tree tree) {
        double length = 0;
        for (Node node : tree.getExternalNodes()) {
                length += node.getLength();
        }
        return length;
    }

    private double getBranchLengthSD(Tree tree, double length) {
    	double variance = 0;
        double mean = length / (tree.getInternalNodeCount() - 1);
    	for (Node node : tree.getInternalNodes()) {
    		if (!node.isRoot()) {
    			variance += Math.pow(node.getLength() - mean, 2);
    		}
    	}
        variance /= (tree.getInternalNodeCount() - 1);
        return Math.sqrt(variance);
    }

	@Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        return 3; // height, length, variance
    }

    @Override
    public double getArrayValue() {
        return treeInput.get().getRoot().getHeight();
    }

    @Override
    public double getArrayValue(int dim) {
    	if (dim == 0) {
    		return treeInput.get().getRoot().getHeight();
    	}
    	return getLength(treeInput.get());
    }
}
