package adb.distribution;

import adb.util.Utils;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;


@Description("Likelihood of Tree under ADB model.")
// TODO: figure out tree factors for multi-type and annotated trees...
public class ADBTreeDistribution extends SpeciesTreeDistribution {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);

    // computing options
    public Input<Integer> maxIterationsInput =
            new Input<>("maxIterations", "maximum number of iterations for numerical integration",100);
    public Input<Double> toleranceInput =
            new Input<>("tolerance", "tolerance for numerical integration",1e-12);
    public Input<Integer> nStepsInput =
            new Input<>("nSteps", "number of time steps for FFT", (int)Math.pow(2, 14));
    public Input<Boolean> approxInput =
            new Input<>("approx", "approximate branch probabilities (default true)", true);
    public Input<Boolean> conditionOnOriginInput =
            new Input<>("conditionOnOrigin", "condition on time since origin otherwise on root height (default true)", true);
    public Input<Boolean> useAnalyticalBDSolutionInput =
            new Input<>("useAnalyticalBDSolution", "use analytical solution if shape is 1 (default false)", false);


    private Parameterization parameterization;
    private TreeInterface tree;

    // options
    private int maxIt;
    private double tol;
    private int nSteps;
    private boolean approx;
    private boolean conditionOnOrigin;
    private boolean useBD;


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        parameterization = parameterizationInput.get();
        tree = treeInput.get();

        // make sure that all tips are at the same height
        if (!TreeUtils.isUltrametric(tree)) {
            Log.warning("WARNING: This model (tree prior) cannot handle dated tips.");
        }

        maxIt = maxIterationsInput.get();
        tol = toleranceInput.get();
        nSteps = nStepsInput.get();
        approx = approxInput.get();
        conditionOnOrigin = conditionOnOriginInput.get();
        useBD = useAnalyticalBDSolutionInput.get();

        if (approx) {
            if (parameterization.nTypes > 1) {
                throw new IllegalArgumentException("The approximation only works for a single-type branching process!");
            }
            // check whether the shape parameter is an integer
            if (parameterization.shape[0] % 1 != 0) { // alternatively, accept some deviation and round
                throw new IllegalArgumentException("The approximation only works for an integer shape parameter!");
            };
        }

        if (useBD && parameterization.nTypes > 1) {
            throw new IllegalArgumentException("This model only uses the BD analytical solution for a single-type branching process.");
        }

        // check that origin is given if conditioning on it (and not root)
        if (conditionOnOrigin && parameterization.originTime == null) {
            throw new IllegalArgumentException("Please provide the time of origin of the process, or set conditionOnOrigin to false.");
        }

        // check that step size is a power of 2 (required for FFT)
        if (!Utils.isPowerOfTwo(nSteps)) {
            throw new IllegalArgumentException("stepSize must be a power of 2!");
        }

    }


    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        return 0;
    }
}
