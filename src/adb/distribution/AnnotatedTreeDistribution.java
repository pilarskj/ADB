package adb.distribution;

import adb.tree.AnnotatedTree;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;

import java.util.List;
import java.util.Random;

public abstract class AnnotatedTreeDistribution extends Distribution {

    public Input<AnnotatedTree> annotatedTreeInput = new Input<>(
            "annotatedTree", "annotated tree", Input.Validate.REQUIRED);

    protected AnnotatedTree annotatedTree;

    @Override
    public void initAndValidate() {
        annotatedTree = annotatedTreeInput.get();
    }

    // Interface requirements:
    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) { }

}