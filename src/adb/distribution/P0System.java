package adb.distribution;

import adb.distribution.LifetimeDistributions.*;
import adb.util.Utils;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;

import java.util.stream.IntStream;


public class P0System extends CalculationNode {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);

    public Input<LifetimeDistributions> lifetimeDistributionsInput =
            new Input<>("lifetimeDistributions", "", Input.Validate.REQUIRED);

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

    boolean dirty;

    int maxIt;
    double tol;
    int nSteps;
    double[] timeArray;
    double timeStep;

    // stores P0 for each type
    double[][] P0;
    double[][] storedP0;


    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        lifetimeDistributions = lifetimeDistributionsInput.get();

        maxIt = maxIterationsInput.get();
        tol = toleranceInput.get();
        timeArray = timeArrayInput.get();
        timeStep = timeStepInput.get();
        nSteps = timeArray.length;

        P0 = new double[nTypes][nSteps];
        storedP0 = new double[nTypes][nSteps];
        dirty = true;
        calcP0();
    }


    private void calcP0() {

        if (!dirty) return;

        LifetimeDistribution[] distributions = lifetimeDistributions.getLifetimeDistributions();
        double[] d = new double[nTypes];
        double[] rho = new double[nTypes];
        for (int i = 0; i < nTypes; i++) {
            d[i] = parameterization.getDeath(i);
            rho[i] = parameterization.getSampling(i);
        }

        // notation: it = iteration, w = integration variable (time), i,j,k = types
        // initialize matrix
        double[][] X0 = new double[nTypes][nSteps];
        IntStream.range(0, nTypes)
                .parallel()
                .forEach(i -> {
                    double[] cdf = distributions[i].getCDF();
                    for (int w = 0; w < nSteps; w++) {
                        X0[i][w] = (1 - rho[i]) * (1 - cdf[w]) + d[i] * cdf[w];
                    }
                });

        // set up iteration
        double err = 1;
        int it = 0;
        double[][] X = X0;

        // iterate
        while (err > tol && it < maxIt) {
            double[][] Xi = new double[nTypes][nSteps];

            for (int i = 0; i < nTypes; i++) {
                // get vectors for convolution
                double[] y = new double[nSteps];
                for (int w = 0; w < nSteps; w++) { // multiply elementwise on times
                    for (int j = 0; j < nTypes; j++) { // sum over all types k
                        y[w] += parameterization.getSymTransition(i,j) * X[j][w] * X[j][w] +
                                parameterization.getAsymTransition(i,j) * X[i][w] * X[j][w];
                    }
                }
                // partially convolve
                double[] I = Utils.convolveFFT(distributions[i].getTransformedPDF(), y, nSteps, timeStep);
                // sum
                for (int w = 0; w < nSteps; w++) {
                    Xi[i][w] = X0[i][w] + (1 - d[i]) * I[w];
                }
            }

            // compute error
            err = Utils.getError(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxIt) {
            System.err.printf("calculateP0 Warning: max iterations reached with error: %.2f%n", err);
        }

        // set final
        P0 = X;
        dirty = false;
    }


    public double[][] getP0() {
        calcP0();
        return P0;
    }


    @Override
    public boolean requiresRecalculation() {
        dirty = true;
        return true;
    }


    @Override
    protected void store() {
        for (int i = 0; i < nTypes; i++) {
            System.arraycopy(P0[i], 0, storedP0[i], 0, nSteps);
        }
        super.store();
    }


    @Override
    protected void restore() {
        double[][] tmp;
        tmp = P0;
        P0 = storedP0;
        storedP0 = tmp;
        super.restore();
    }

}
