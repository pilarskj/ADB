package adb.distribution;

import adb.tree.AnnotatedNode;
import adb.tree.AnnotatedTree;
import adb.util.Utils;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.util.Pair;

import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

import static adb.util.Utils.TRANSFORM_FORWARD;
import static adb.tree.AnnotatedNode.getType;


@Description("Likelihood of a tree under ADB model.")
// TODO: assert tree factors for multi-type and annotated trees
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
    public Input<Boolean> conditionOnOriginInput =
            new Input<>("conditionOnOrigin", "condition on time since origin otherwise on root height (default true)", true);
    public Input<Boolean> approxInput =
            new Input<>("approx", "approximate branch probabilities (default true)", true);
    public Input<Boolean> useAnalyticalBDSolutionInput =
            new Input<>("useAnalyticalBDSolution", "use analytical solution if shape is 1 (default false)", false);


    Parameterization parameterization;
    TreeInterface tree;

    // options
    int maxIt;
    double tol;
    int nSteps;
    boolean conditionOnOrigin;
    boolean approx;
    boolean useBD;

    // time array
    double[] timeArray;
    double timeStep;

    int nTypes;

    // calculation nodes
    LifetimeDistributions lifetimeDistributions;
    P0System P0System;
    HashMap<Integer, UnivariateFunction> P0Map;

    // optional
    P1System P1System;
    HashMap<Pair<Integer,Integer>, UnivariateFunction> P1Map;

    // for approximation:
    // use cashing to avoid creating a new GammaDistribution object with shape i*b for i=1,2,... in each iteration
    // create a map to store the distributions for each i (once needed) and re-use them
    // use a thread-safe and dynamic map
    ConcurrentHashMap<Integer, GammaDistribution> gammaCache = new ConcurrentHashMap<>();


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        tree = treeInput.get();

        // make sure that all tips are at the same height
        if (!TreeUtils.isUltrametric(tree)) {
            Log.warning("WARNING: This model (tree prior) cannot handle dated tips.");
        }

        // make sure that all tips have valid types
        for (Node tip : tree.getExternalNodes()) {
            int type = getType(tip);
            if (type < 0 || type >= nTypes) {
                throw new IllegalArgumentException("Not all types at tips are valid, please check the input.");
            }
        }

        maxIt = maxIterationsInput.get();
        tol = toleranceInput.get();
        nSteps = nStepsInput.get();
        conditionOnOrigin = conditionOnOriginInput.get();
        approx = approxInput.get();
        useBD = useAnalyticalBDSolutionInput.get();

        if (approx) {
            if (nTypes > 1) {
                throw new IllegalArgumentException("The approximation only works for a single-type branching process!");
            }
            // check whether the shape parameter is an integer
            if (!parameterization.shapeIsInteger()) {
                throw new IllegalArgumentException("The approximation only works for an integer shape parameter!");
            };
        }

        if (useBD && nTypes > 1) {
            throw new IllegalArgumentException("This model only uses the BD analytical solution for a single-type branching process.");
        }

        // check that step size is a power of 2 (required for FFT)
        if (!Utils.isPowerOfTwo(nSteps)) {
            throw new IllegalArgumentException("stepSize must be a power of 2!");
        }

        // create time array for computation
        if (conditionOnOrigin) {
            if (parameterization.getOriginTime() == null) {
                throw new IllegalArgumentException("Please provide the time of origin of the process, or set conditionOnOrigin to false.");
            } else {
                timeStep = parameterization.getOriginTime() / nSteps;
                timeArray = Utils.linSpace(0, parameterization.getOriginTime(), nSteps);
            }
        } else {
            timeStep = tree.getRoot().getHeight() / nSteps;
            timeArray = Utils.linSpace(0, tree.getRoot().getHeight(), nSteps);
        }

        // initialize internal calculation nodes
        lifetimeDistributions = new LifetimeDistributions();
        lifetimeDistributions.init(
                parameterization.getLifetimeParameter(),
                parameterization.getShapeParameter(),
                timeArray
        );

        P0System = new P0System();
        P0System.init(
                parameterization,
                lifetimeDistributions,
                maxIt, tol, timeArray, timeStep
        );

        if (!(tree instanceof AnnotatedTree)) { // TODO: or xml input tag?
            P1System.init(
                    parameterization,
                    lifetimeDistributions,
                    P0System,
                    maxIt, tol, timeArray, timeStep
            );
        }

    }


    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        double[][] P0 = P0System.getP0();
        // stop if extinction is certain
        if (P0[parameterization.getOriginType()][nSteps - 1] == 1.0) {
            return Double.NEGATIVE_INFINITY;
        }

        // extend time array to calculate probabilities of tiny branches
        double[] extTimeArray = new double[nSteps + 1];
        extTimeArray[0] = 0;
        System.arraycopy(timeArray, 0, extTimeArray, 1, nSteps);

        // extend and interpolate P0
        P0Map = new HashMap<>();
        for (int i = 0; i < nTypes; i++) {
            double[] extP0 = new double[nSteps + 1];
            extP0[0] = 1 - parameterization.getSampling(i);
            System.arraycopy(extP0, 0, P0[i], 1, nSteps);
            UnivariateFunction function = Utils.interpolator.interpolate(extTimeArray, extP0);
            P0Map.put(i, function);
        }

        if (!(tree instanceof AnnotatedTree)) {
            double[][][] P1 = P1System.getP1();

            // extend and interpolate P1
            P1Map = new HashMap<>();
            for (int i = 0; i < nTypes; i++) {
                for (int j = 0; j < nTypes; j++) {
                    double[] extP1 = new double[nSteps + 1];
                    if (i == j) { extP1[0] = parameterization.getSampling(i); } else { extP1[0] = 0; }
                    System.arraycopy(extP1, 0, P1[i][j], 1, nSteps);
                    UnivariateFunction function = Utils.interpolator.interpolate(extTimeArray, extP1);
                    P1Map.put(new Pair<>(i, j), function);
                }
            }
        }

        // TODO: Tree likelihood comes here!

        return 0;
    }


    // track time backwards (start < end)
    private double[] calculateSubtreeLikelihood(Node node, double start, double end) {

        double[] likelihood = new double[nTypes];
        int type = getType(node);

        // upstream branch
        double[][] branchDensity = calculateNodeLikelihood(node, start, end);

        if (node.isLeaf()) {
            for (int i = 0; i < nTypes; i++) {
                likelihood[i] = branchDensity[i][type];
            }

        } else {
            Node leftChild = node.getLeft();
            Node rightChild = node.getRight();
            double[] leftSubtreeLik = calculateSubtreeLikelihood(leftChild, leftChild.getHeight(), start);
            double[] rightSubtreeLik = calculateSubtreeLikelihood(rightChild, rightChild.getHeight(), start);

            for (int i = 0; i < nTypes; i++) {
                double lik = 0;
                for (int j = 0; j < nTypes; j++) {
                    double sum = 0;
                    for (int k = 0; k < nTypes; k++) {
                        sum += parameterization.getSymTransition(j, k) * leftSubtreeLik[k] * rightSubtreeLik[k] +
                                0.5 * parameterization.getAsymTransition(j, k) *
                                        (leftSubtreeLik[j] * rightSubtreeLik[k] + leftSubtreeLik[k] * rightSubtreeLik[j]);
                    }
                    lik += branchDensity[i][j] * sum;
                }
                likelihood[i] = lik;
            }
        }

        return likelihood;
    }


    // calculate the density of the branch upstream of node
    private double[][] calculateNodeLikelihood(Node node, double start, double end) {

        double[][] likelihood = new double[nTypes][nTypes];

        if (node instanceof AnnotatedNode) {
            return likelihood; // TODO: branch segments

        } else {
            if (node.isLeaf()) { // tip - sampling event
                assert start == 0.0;
                for (int i = 0; i < nTypes; i++) {
                    for (int j = 0; j < nTypes; j++) {
                        likelihood[i][j] = P1Map.get(new Pair<>(i, j)).value(end);
                    }
                }

            } else {
                if (approx) {
                    likelihood[0][0] = approximateBranchDensity(start, end);
                } else {
                    likelihood = calculateBranchDensity(start, end);
                }
            }

            return likelihood;
        }
    }


    // calculate the density of an internal branch (s,e) with types i -> j
    // requires solving a system of integral equations; very slow
    private double[][] calculateBranchDensity(double start, double end) {

        double[][] B = new double[nTypes][nTypes];

        // generate linearly spaced values between start and end
        double[] seq = new double[nSteps];
        double[] ageSeq = new double[nSteps];
        double dx = (end - start) / nSteps;
        for (int w = 0; w < nSteps; w++) {
            seq[w] = start + dx * (w + 1);
            ageSeq[w] = seq[w] - start;
        }

        // calculate the FFT from PDF of the gamma distribution and interpolate P0 per type, initialize matrix
        Complex[][] pdfFFT = new Complex[nTypes][nSteps * 2];
        double[][] P0 = new double[nTypes][nSteps];
        double[][][] X0 = new double[nTypes][nTypes][nSteps];
        IntStream.range(0, nTypes)
                .parallel()
                .forEach(i -> {
                    double[] pdf = new double[nSteps];
                    double lifetime = parameterization.getLifetime(i);
                    double shape = parameterization.getShape(i);
                    double scale = lifetime / shape;
                    GammaDistribution gammaDist = new GammaDistribution(shape, scale);
                    for (int w = 0; w < nSteps; w++) {
                        pdf[w] = Math.exp(gammaDist.logDensity(ageSeq[w])); // get density
                        P0[i][w] = P0Map.get(i).value(seq[w]); // extrapolate P0
                        X0[i][i][w] = (1 - parameterization.getDeath(i)) * pdf[w]; // initialize matrix
                    }

                    Complex[] Ft = Utils.fft.transform(Utils.padZeros(pdf), TRANSFORM_FORWARD);
                    for (int w = 0; w < nSteps * 2; w++) {
                        pdfFFT[i][w] = Ft[w];
                    }
                });

        // set up iteration
        double err = 1;
        int it = 0;
        double[][][] X = X0;

        // iterate
        while (err > tol && it < maxIt) {
            double[][][] Xi = new double[nTypes][nTypes][nSteps];

            for (int i = 0; i < nTypes; i++) {
                for (int j = 0; j < nTypes; j++) {
                    // get vectors for convolution
                    double[] y = new double[nSteps];
                    for (int w = 0; w < nSteps; w++) { // multiply elementwise on times
                        for (int k = 0; k < nTypes; k++) { // sum over all types k
                            y[w] += parameterization.getSymTransition(i, k) * P0[k][w] * X[k][j][w] +
                                    0.5 * parameterization.getAsymTransition(i, k) * (P0[i][w] * X[k][j][w] + P0[k][w] * X[i][j][w]);
                        }
                    }

                    // extract column from the pdf matrix
                    Complex[] Ft = pdfFFT[i];

                    // partially convolve
                    double[] I = Utils.convolveFFT(Ft, y, nSteps, dx);

                    // sum
                    for (int w = 0; w < nSteps; w++) {
                        Xi[w][i][j] = X0[w][i][j] + 2 * (1 - parameterization.getDeath(i)) * I[w];
                    }
                }
            }

            // compute error
            err = Utils.getMatrixError3D(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxIt) {
            System.err.printf("calculateBranchDensity Warning: max iterations reached with error: %.2f%n", err);
        }

        // take final time slice only
        for (int i = 0; i < nTypes; i++) {
            for (int j = 0; j < nTypes; j++) {
                B[i][j] = X[i][j][nSteps - 1];
            }
        }

        return B;
    }


    // single-type approximation
    private double approximateBranchDensity(double start, double end) {

        assert nTypes == 1;
        double[] P0 = P0System.getP0()[0];

        // get average P0 over branch
        // use binary search to find closest indices to the branch lengths in time array (t0 is sorted per definition!)
        double[] P0Slice = Arrays.copyOfRange(P0, Utils.findClosestIndex(timeArray, start), Utils.findClosestIndex(timeArray, end) + 1);
        double P0M = Utils.getMean(P0Slice); // TODO: take extrapolated values using P0Map

        // initialize the approximation
        double lifetime = parameterization.getLifetime(0);
        double shape = parameterization.getShape(0);
        double scale = lifetime / shape;
        int k = (int) ((end - start) / (shape * scale)); // around this k, the term b_k will be maximal // TODO: update maximum (factor)
        double density = 0;

        int i = 0; // increase (start with 0 and treat k as minimum number of terms)
        double term = 1;
        while (term > tol || i <= k) { // only stops when both term <= tol and i > k!
            GammaDistribution gammaDist;
            // check if GammaDistribution with the given shape is already in the cache
            if (gammaCache.containsKey(i + 1)) {
                gammaDist = gammaCache.get(i + 1);
            } else { // otherwise, add it
                gammaDist = new GammaDistribution((i + 1) * shape, scale);
                gammaCache.put(i + 1, gammaDist);
            }
            term = Math.pow(2, i) * Math.pow(1 - parameterization.getDeath(0), i + 1) * Math.pow(P0M, i) * Math.exp(gammaDist.logDensity(end - start));
            density += term;
            i++;
        }

        return density;
    }

}
