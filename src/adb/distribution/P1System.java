package adb.distribution;

import adb.distribution.LifetimeDistributions.LifetimeDistribution;
import adb.util.Utils;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;

import java.util.stream.IntStream;


public class P1System extends CalculationNode {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);

    public Input<LifetimeDistributions> lifetimeDistributionsInput =
            new Input<>("lifetimeDistributions", "", Input.Validate.REQUIRED);

    public Input<P0System> P0SystemInput =
            new Input<>("P0System", "", Input.Validate.REQUIRED);

    public Input<Integer> maxIterationsInput =
            new Input<>("maxIterations", "",Input.Validate.REQUIRED);
    public Input<Double> toleranceInput =
            new Input<>("tolerance", "", Input.Validate.REQUIRED);
    public Input<double[]> timeArrayInput =
            new Input<>("timeArray", "", Input.Validate.REQUIRED);
    public Input<Double> timeStepInput =
            new Input<>("timeStep", "", Input.Validate.REQUIRED);

    Parameterization parameterization;
    int nTypes;
    LifetimeDistributions lifetimeDistributions;
    P0System P0System;

    boolean dirty;

    int maxIt;
    double tol;
    int nSteps;
    double[] timeArray;
    double timeStep;

    // stores P1 for each type
    double[][][] P1;
    double[][][] storedP1;


    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        lifetimeDistributions = lifetimeDistributionsInput.get();
        P0System = P0SystemInput.get();

        maxIt = maxIterationsInput.get();
        tol = toleranceInput.get();
        timeArray = timeArrayInput.get();
        timeStep = timeStepInput.get();
        nSteps = timeArray.length;

        P1 = new double[nTypes][nTypes][nSteps];
        storedP1 = new double[nTypes][nTypes][nSteps];
        dirty = true;
        calculateP1();
    }

    private void calculateP1() {

        if (!dirty) return;

        LifetimeDistribution[] distributions = lifetimeDistributions.getLifetimeDistributions();
        double[] d = new double[nTypes];
        double[] rho = new double[nTypes];
        for (int i = 0; i < nTypes; i++) {
            d[i] = parameterization.getDeath(i);
            rho[i] = parameterization.getSampling(i);
        }

        double[][] P0 = P0System.getP0();


        // notation: it = iteration, w = integration variable (time), i,j,k = types
        // initialize matrix
        double[][][] X0 = new double[nTypes][nTypes][nSteps];
        IntStream.range(0, nTypes)
                .parallel()
                .forEach(i -> {
                    double[] cdf =  distributions[i].getCDF();
                    for (int w = 0; w < nSteps; w++) {
                        X0[i][i][w] = rho[i] * (1 - cdf[w]);
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
                                    0.5 * parameterization.getSymTransition(i, k) * (P0[i][w] * X[k][j][w] + P0[k][w] * X[i][j][w]);
                        }
                    }

                    // partially convolve
                    double[] I = Utils.convolveFFT(distributions[i].getTransformedPDF(), y, nSteps, timeStep);

                    // sum
                    for (int w = 0; w < nSteps; w++) {
                        Xi[i][j][w] = X0[i][j][w] + 2 * (1 - d[i]) * I[w];
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
            System.err.printf("calculateP1 Warning: max iterations reached with error: %.2f%n", err);
        }

        // set state
        P1 = X;

        dirty = false;
    }


    public double[][][] getP1() {
        calculateP1();
        return P1;
    }


    @Override
    public boolean requiresRecalculation() {
        dirty = true;
        return true;
    }


    @Override
    protected void store() {
        for (int i = 0; i < nTypes; i++) {
            for (int j = 0; j < nTypes; i++) {
                System.arraycopy(P1[i][j], 0, storedP1[i][j], 0, nSteps);
            }
        }
        super.store();
    }


    @Override
    protected void restore() {
        double[][][] tmp;
        tmp = P1;
        P1 = storedP1;
        storedP1 = tmp;
        super.restore();
    }


}
