package adb.util;

import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedTree;
import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.io.PrintStream;

public class BranchingTreeLogger extends BEASTObject implements Loggable {

    public Input<AnnotatedTree> annotatedTreeInput =
            new Input<>("annotatedTree", "annotated tree to log as standard branching tree", Input.Validate.REQUIRED);

    public Input<Boolean> printMetaDataInput =
            new Input<>("printMetaData", "print node metadata (including type) (default true)", true);

    Tree tree;
    boolean printMetaData;

    @Override
    public void initAndValidate() {
        tree = annotatedTreeInput.get();
        printMetaData = printMetaDataInput.get();
    }

    @Override
    public void init(PrintStream out) {
        tree.init(out);
    }

    @Override
    public void log(long nSample, PrintStream out) {

        if (printMetaData) {
            // set up metadata string
            for (Node node : tree.getNodesAsArray()) {
                node.setMetaData("type", ((AnnotatedNode)node).getType());
                node.metaDataString = String.format("%s=%d", "type", ((AnnotatedNode)node).getType());
            }
        }
        out.print("tree STATE_" + nSample + " = ");
        out.print(tree.getRoot().toSortedNewick(new int[1], printMetaData));
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        tree.close(out);
    }

}
