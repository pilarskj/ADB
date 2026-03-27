package adb.distribution;

import adb.util.Utils;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.complex.Complex;

import java.util.HashMap;
import java.util.List;
import java.util.stream.IntStream;

public class P0System extends CalculationNode {

    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "ADB parameterization", Input.Validate.REQUIRED);

    Parameterization parameterization;
    List<LifetimeDistribution> lifetimeDist;
    int nTypes;

    private double maxIt;
    private double tol;
    int nSteps;
    private double[] seq;
    private double dx;


    // stores P0 for each type as a function of time
    public HashMap<Integer, UnivariateFunction> P0Map;


    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        // get number of types and time steps
        nTypes = parameterization.nTypes;
    }

    // TODO: adapt function for calculating the extinction probability
    // or return map?
    /* public void calcP0() {
        Complex[][] pdfFFT, double[][] cdf, double[] seq, double dx,
        double[] d, double rho, double[][] Xsi_s, double[][] Xsi_as,
        int maxIt, double tol

        // notation: it = iteration, w = integration variable (time), i,j,k = types
        // initialize matrix
        double[][] X0 = new double[nSteps][nTypes];
        IntStream.range(0, nTypes)
                .parallel()
                .forEach(i -> {
                    for (int w = 0; w < nSteps; w++) {
                        X0[w][i] = (1 - rho) * (1 - cdf[w][i]) + d[i] * cdf[w][i];
                    }
                });

        // set up iteration
        double err = 1;
        int it = 0;
        double[][] X = X0;

        // iterate
        while (err > tol && it < maxIt) {
            double[][] Xi = new double[nSteps][nTypes];

            for (int i = 0; i < nTypes; i++) {
                // get vectors for convolution
                double[] y = new double[nSteps];
                for (int w = 0; w < nSteps; w++) { // multiply elementwise on times
                    for (int j = 0; j < nTypes; j++) { // sum over all types k
                        y[w] += Xsi_s[i][j] * X[w][j] * X[w][j] + Xsi_as[i][j] * X[w][i] * X[w][j];
                    }
                }
                // extract column from the pdf matrix
                Complex[] Ft = new Complex[nSteps*2];
                for (int w = 0; w < nSteps*2; w++) {
                    Ft[w] = pdfFFT[w][i];
                }

                // partially convolve
                double[] I = Utils.convolveFFT(Ft, y, nSteps, dx);

                // sum
                for (int w = 0; w < nSteps; w++) {
                    Xi[w][i] = X0[w][i] + (1 - d[i]) * I[w];
                }
            }

            // compute error
            err = Utils.getMatrixError(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxIt) {
            System.err.printf("calcP0 Warning: max iterations reached with error: %.2f%n", err);
        }

        return X;
    } */

}
