package adb;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.Pair;

import java.util.HashMap;
import java.util.stream.IntStream;

import static adb.GammaLogLikelihood.convolveFFT;
import static adb.GammaLogLikelihood.padZeros;
import static adb.MTLogLikelihood.calcMTB;
import static adb.MTLogLikelihood.calcMTP0;

/*
Class for solving the equations for P0, P1 and B in the multi-type case
and for calculating the likelihood of an uncoloured (or coloured) tree
based on the parameters, branching times, and types at tips (and internal nodes).
 */
public class MTLogLikelihood3D {

    public static FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
    public static LinearInterpolator interpolator = new LinearInterpolator();

    public static double calcMTLogLikelihood3D(double[] a, double[] b, double[] d, double rho,
                                               double[][] Xsi_as, double[][] Xsi_s,
                                               double t_or, int type_or,
                                               BranchList branches,
                                               int maxIt, double tolP, double tolB, int mP, int mB) {

        // get number of types
        int ntypes = a.length;
        assert b.length == ntypes;
        assert d.length == ntypes;

        // assert type at origin
        branches.getBranchByIndex(0).endType = type_or;

        // generate linearly spaced values between 0 and origin
        int m = mP; // the number of time steps must be a power of 2 (required by FFT!)
        double[] tSeq = new double[m];
        double dx = t_or / m;
        for (int i = 0; i < m; i++) {
            tSeq[i] = dx * (i + 1);
        }
        assert tSeq[m - 1] == t_or;

        // calculate CDF and FFT from PDF from 0 to origin
        double[][] cdf = new double[m][ntypes];
        Complex[][] pdfFFT = new Complex[m*2][ntypes];
        IntStream.range(0, ntypes)
                .parallel()
                .forEach(x -> {
                    // get density and cumulative probability
                    double[] pdf = new double[m];
                    GammaDistribution gammaDist = new GammaDistribution(b[x], a[x]);
                    for (int i = 0; i < m; i++) {
                        pdf[i] = Math.exp(gammaDist.logDensity(tSeq[i]));
                        cdf[i][x] = gammaDist.cumulativeProbability(tSeq[i]);
                    }

                    // perform FFT
                    Complex[] Ft = fft.transform(padZeros(pdf), TransformType.FORWARD);
                    for (int i = 0; i < m*2; i++) {
                        pdfFFT[i][x] = Ft[i];
                    }
                });

        // calculate extinction probability over time
        double[][] P0 = calcMTP0(pdfFFT, cdf, d, rho, Xsi_as, Xsi_s, tSeq, dx, maxIt, tolP);

        // calculate probability of single descendant over time
        double[][][] P1 = calcMTP13D(pdfFFT, cdf, d, rho, Xsi_as, Xsi_s, tSeq, P0, dx, maxIt, tolP);

        // extend tSeq to calculate probabilities of tiny branches
        double[] extSeq = new double[tSeq.length + 1];
        extSeq[0] = 0;
        System.arraycopy(tSeq, 0, extSeq, 1, tSeq.length);

        // interpolate P1
        HashMap<Pair<Integer,Integer>, UnivariateFunction> P1Map = new HashMap<>();
        for (int i = 0; i < ntypes; i++) {
            for (int j = 0; j < ntypes; j++) {
                // extract array from P1 (and extend to 0)
                double[] X = new double[m + 1];
                if (i == j) { X[0] = rho; } else { X[0] = 0; }
                for (int w = 0; w < m ; w++) {
                    X[w + 1] = P1[w][i][j];
                }
                // interpolate
                UnivariateFunction function = interpolator.interpolate(extSeq, X);
                P1Map.put(new Pair<>(i, j), function);
            }
        }

        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();

        // collect times of internal branches and relationships for recursion
        double[] intS = new double[n_int];
        double[] intE = new double[n_int];
        int[] left_child = new int[n_int];
        int[] right_child = new int[n_int];
        for (int i = 0; i < n_int; i++) {
            Branch branch = branches.getBranchByIndex(i);
            intS[i] = branch.startTime;
            intE[i] = branch.endTime;
            left_child[i] = branch.leftIndex;
            right_child[i] = branch.rightIndex;
        }

        /*// sum probabilities over external branches
        double logP1 = 0;
        for (int x = 0; x < n_ext; x++) {
            Branch branch = branches.getBranchByIndex(x + n_int);
            double p = 0;
            for (int i = 0; i < ntypes; i++) {
                p += P1Map.get(new Pair<>(i, branch.startType)).value(branch.endTime);
            }
            // double p = functionMap.get(new Pair<>(branch.endType, branch.startType)).value(branch.endTime);
            logP1 += Math.log(p); // if probability at some tip < 0 due to rounding errors, this will be -Infinity
        }
        if (logP1 == Double.NEGATIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }
        return logP1; */

        // calculate probabilities of internal branches
        double[][] B = calcMTB(a, b, d, Xsi_as, Xsi_s, intS, intE, tSeq, P0, maxIt, tolB, mB);

        // start recursion
        double treeL = calcLikelihoodTipTyped(1, branches, n_ext, ntypes, 0, type_or,
                left_child, right_child, Xsi_as, Xsi_s, P1Map, B);

        return Math.log(treeL);
    }

    // Recursive function for calculating the tree likelihood (for tip-typed tree)
    public static double calcLikelihoodTipTyped(double likelihood, BranchList branches, int ntips, int ntypes,
                                                int stem_index, int stem_type, int[] left_child, int[] right_child,
                                                double[][] Xsi_as, double[][] Xsi_s,
                                                HashMap<Pair<Integer,Integer>, UnivariateFunction> P1Map, double[][] B) {

        // at tips
        if (stem_index >= ntips - 1) {
            Branch branch = branches.getBranchByIndex(stem_index);
            return P1Map.get(new Pair<>(stem_type, branch.startType)).value(branch.endTime);

        } else {
            // recursion
            double subtreeL = 0;
            // loop over all possible types at internal nodes
            for (int j = 0; j < ntypes; j++){

                // get subtree (new stems)
                int right_stem = right_child[stem_index];
                int left_stem = left_child[stem_index];

                subtreeL += Xsi_s[stem_type][j] *
                        calcLikelihoodTipTyped(likelihood, branches, ntips, ntypes, left_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1Map, B) *
                        calcLikelihoodTipTyped(likelihood, branches, ntips, ntypes, right_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1Map, B) +
                        Xsi_as[stem_type][j] *
                                (calcLikelihoodTipTyped(likelihood, branches, ntips, ntypes, left_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1Map, B) *
                                        calcLikelihoodTipTyped(likelihood, branches, ntips, ntypes, right_stem, stem_type, left_child, right_child, Xsi_as, Xsi_s, P1Map, B) +
                                        calcLikelihoodTipTyped(likelihood, branches, ntips, ntypes, left_stem, stem_type, left_child, right_child, Xsi_as, Xsi_s, P1Map, B) *
                                                calcLikelihoodTipTyped(likelihood, branches, ntips, ntypes, right_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1Map, B));
                        ;
            }

            likelihood = B[stem_index][stem_type] * subtreeL;
            return likelihood;
        }
    }


    // Function for calculating the probability of a single descendant (considering start and end type)
    public static double[][][] calcMTP13D(Complex[][] pdfFFT, double[][] cdf, double[] d, double rho, double[][] Xsi_as, double[][] Xsi_s,
                                          double[] t0, double[][] P0, double dx, int maxIt, double tol) {

        // get number of types and time steps
        int ntypes = cdf[0].length;
        int m = t0.length;

        // initialize matrix
        double[][][] X0 = new double[m][ntypes][ntypes];
        for (int w = 0; w < m; w++) {
            for (int i = 0; i < ntypes; i++) {
                X0[w][i][i] = rho * (1 - cdf[w][i]);
            }
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[][][] X = X0;

        // iterate
        while (err > tol && it < maxIt) {
            double[][][] Xi = new double[m][ntypes][ntypes];

            for (int i = 0; i < ntypes; i++) {
                for (int j = 0; j < ntypes; j++) {
                    // get vectors for convolution
                    double[] y = new double[m];
                    for (int w = 0; w < m; w++) { // multiply elementwise on times
                        for (int k = 0; k < ntypes; k++) { // sum over all types k
                            y[w] += Xsi_s[i][k] * P0[w][k] * X[w][k][j] +
                                    Xsi_as[i][k] * (P0[w][i] * X[w][k][j] + P0[w][k] * X[w][i][j]);
                        }
                    }

                    // extract column from the pdf matrix
                    Complex[] Ft = new Complex[m * 2];
                    for (int w = 0; w < m * 2; w++) {
                        Ft[w] = pdfFFT[w][i];
                    }

                    // partially convolve
                    double[] I = convolveFFT(Ft, y, m, dx);

                    // sum
                    for (int w = 0; w < m; w++) {
                        Xi[w][i][j] = X0[w][i][j] + 2 * (1 - d[i]) * I[w];
                    }
                }

                // compute error
                err = getMatrixError3D(X, Xi);

                // update
                X = Xi;
                it++;
            }

            if (it == maxIt) {
                System.err.printf("calcP1 Warning: max iterations reached with error: %.2f%n", err);
            }
        }

        return X;
    }


    // Helper method for calculating the distance between two 3D arrays
    private static double getMatrixError3D(double[][][] X, double[][][] Y){
        int n = X.length; // number of arrays
        int m = X[0].length;  // number of rows
        int o = X[0][0].length;  // number of columns

        double maxSum = 0;

        // iterate over arrays
        for (int i = 0; i < n; i++) {
            double sum = 0;

            // substract matrices element-wise
            for (int row = 0; row < m; row++) {
                for (int col = 0; col < o; col++) {
                    // sum the absolute differences for each element
                    sum += Math.abs(X[i][row][col] - Y[i][row][col]);
                }
            }

            // update maxSum if this array has a larger sum
            if (sum > maxSum) {
                maxSum = sum;
            }
        }

        return maxSum;
    }

}



