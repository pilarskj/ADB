package adb.distribution;

import adb.tree.AnnotatedTree;
import beast.base.core.Description;
import beast.base.core.Input;

@Description("Likelihood of AnnotatedTree under ADB model.")
public class ADBAnnotatedTreeDistribution extends AnnotatedTreeDistribution {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);
    // + all options as Inputs

    protected Parameterization parameterization;
    protected AnnotatedTree annotatedTree;
    // + all options as variables

    public ADBAnnotatedTreeDistribution() { };


}
