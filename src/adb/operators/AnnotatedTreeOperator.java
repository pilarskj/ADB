package adb.operators;

import adb.tree.AnnotatedTree;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;

@Description("This operator generates proposals for an AnnotatedTree.")
public abstract class AnnotatedTreeOperator extends Operator {

    public Input<AnnotatedTree> annotatedTreeInput = new Input<>(
            "annotatedTree", "Annotated tree on which to operate.",
            Input.Validate.REQUIRED);

    protected AnnotatedTree tree;

    @Override
    public void initAndValidate() {
        tree = annotatedTreeInput.get();
    }

    /* potentially, copy functions from TreeOperator */
}
