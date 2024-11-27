package adbp;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import static adbp.GammaLogLikelihood.convolveFFT;
import static adbp.GammaLogLikelihood.padZeros;

/*
Class for solving the equations for P0, P1 and B in the multi-type case
and for calculating the likelihood of an uncoloured (or coloured) tree
based on the parameters, branching times, and types at tips (and internal nodes).
 */
public class MTLogLikelihood {

    // can those be initialized once and re-used?
    public static FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
    public static LinearInterpolator interpolator = new LinearInterpolator();

    public static double calcMTLogLikelihood(double[] a, double[] b, double[] d, double rho,
                                             double[][] Xsi_as, double[][] Xsi_s,
                                             double t_or, int[] types, //int type_or,
                                             double[] intS, double[] intE, double[] extE,
                                             int[] left_child, int[] right_child,
                                             int maxIt, double tolP, double tolB, int mP, int mB) {

        // get number of types
        int ntypes = a.length;
        assert b.length == ntypes;
        assert d.length == ntypes;

        // get number of tips
        int ntips = extE.length;

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
                        pdf[i] = gammaDist.density(tSeq[i]);
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

        // calculate probability of single descendants at tips
        double[][] P1 = calcMTP1(pdfFFT, cdf, d, rho, Xsi_as, Xsi_s, extE, tSeq, P0, dx, maxIt, tolP);

        // calculate probabilities of internal branches
        double[][] B = calcMTB(a, b, d, Xsi_as, Xsi_s, intS, intE, tSeq, P0, maxIt, tolB, mB);
        // System.out.println(Arrays.deepToString(B));

        // start recursion
        double treeL = calcLikelihoodColoured(1, ntips, 0, types[0], types,
                left_child, right_child, Xsi_as, Xsi_s, P1, B);
        // uncoloured (super slow)
        // double treeL = calcLikelihoodUncoloured(1, ntips, ntypes, 0, type_or, left_child, right_child,
        //        Xsi_as, Xsi_s, P1, B);

        return Math.log(treeL);
    }


    // Recursive function for calculating the tree likelihood (for a coloured tree)
    public static double calcLikelihoodColoured(double likelihood, int ntips, int stem_index, int stem_type,
                                                int[] types, int[] left_child, int[] right_child,
                                                double[][] Xsi_as, double[][] Xsi_s,
                                                double[][] P1, double[][] B) {

        // at tips
        if (stem_index >= ntips - 1) {
            likelihood = P1[stem_index - (ntips - 1)][stem_type];

        } else {
            // get subtree (new stems)
            int right_stem = right_child[stem_index];
            int left_stem = left_child[stem_index];

            int right_type = types[right_stem];
            int left_type = types[right_stem];

            // find transition probability
            double xsi;
            if (right_type == left_type) {
                xsi = Xsi_s[stem_type][right_type];
            } else {
                if (right_type == stem_type) {
                    xsi = 0.5 * Xsi_as[stem_type][left_type];
                } else if (left_type == stem_type) {
                    xsi = 0.5 * Xsi_as[stem_type][right_type];
                } else {
                   xsi = 0;
                }
            }

            // recursion
            double subtreeL = xsi *
                    calcLikelihoodColoured(likelihood, ntips, left_stem, left_type, types, left_child, right_child, Xsi_as, Xsi_s, P1, B) *
                    calcLikelihoodColoured(likelihood, ntips, right_stem, right_type, types, left_child, right_child, Xsi_as, Xsi_s, P1, B);

            likelihood = B[stem_index][stem_type] * subtreeL;
        }

        return likelihood;
    }


    // Recursive function for calculating the tree likelihood (for an uncoloured tree)
    public static double calcLikelihoodUncoloured(double likelihood, int ntips, int ntypes,
                                                  int stem_index, int stem_type, int[] left_child, int[] right_child,
                                                  double[][] Xsi_as, double[][] Xsi_s,
                                                  double[][] P1, double[][] B) {

        // at tips
        if (stem_index >= ntips - 1) {
            return P1[stem_index - (ntips - 1)][stem_type];

        } else {
            // recursion
            double subtreeL = 0;
            // loop over all possible types at internal nodes
            for (int j = 0;  j < ntypes; j++){

                // get subtree (new stems)
                int right_stem = right_child[stem_index];
                int left_stem = left_child[stem_index];

                subtreeL += Xsi_as[stem_type][j] *
                        calcLikelihoodUncoloured(likelihood, ntips, ntypes, right_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1, B) *
                        calcLikelihoodUncoloured(likelihood, ntips, ntypes, left_stem, stem_type, left_child, right_child, Xsi_as, Xsi_s, P1, B) +
                        Xsi_as[stem_type][j] *
                                calcLikelihoodUncoloured(likelihood, ntips, ntypes, right_stem, stem_type, left_child, right_child, Xsi_as, Xsi_s, P1, B) *
                                calcLikelihoodUncoloured(likelihood, ntips, ntypes, left_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1, B) +
                        Xsi_s[stem_type][j] *
                                calcLikelihoodUncoloured(likelihood, ntips, ntypes, right_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1, B) *
                                calcLikelihoodUncoloured(likelihood, ntips, ntypes, left_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1, B);
            }

            likelihood = B[stem_index][stem_type] * subtreeL;
            return likelihood;
        }
    }


    // Function for calculating the extinction probability
    public static double[][] calcMTP0(Complex[][] pdfFFT, double[][] cdf, double[] d, double rho, double[][] Xsi_as, double[][] Xsi_s,
                                      double[] t, double dx, int maxIt, double tol) {

        // get number of types and time steps
        int n = cdf[0].length;
        int m = t.length;

        // initialize matrix
        double[][] X0 = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int x = 0; x < n; x++) {
                X0[i][x] = (1 - rho) * (1 - cdf[i][x]) + d[x] * cdf[i][x];
            }
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[][] X = X0;

        // iterate
        while (err > tol && it < maxIt) {
            double[][] Xi = new double[m][n];

            for (int j = 0; j < n; j++) {
                // get vectors for convolution
                double[] y = new double[m];
                for (int i = 0; i < m; i++) { // multiply elementwise on times
                    for (int k = 0; k < n; k++) { // sum over all types k
                        y[i] += 2 * Xsi_as[j][k] * X[i][j] * X[i][k] + Xsi_s[j][k] * X[i][k] * X[i][k];
                    }
                }
                // extract column from the pdf matrix
                Complex[] Ft = new Complex[m*2];
                for (int i = 0; i < m*2; i++) {
                    Ft[i] = pdfFFT[i][j];
                }

                // partially convolve
                double[] I = convolveFFT(Ft, y, m, dx);

                // sum
                for (int i = 0; i < m; i++) {
                    Xi[i][j] = X0[i][j] + (1 - d[j]) * I[i];
                }
            };

            // compute error
            err = getMatrixError(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxIt) {
            System.err.printf("calcP0 Warning: max iterations reached with error: %.2f%n", err);
        }

        return X;
    }


    // Function for calculating the probability of a single descendant
    public static double[][] calcMTP1(Complex[][] pdfFFT, double[][] cdf, double[] d, double rho, double[][] Xsi_as, double[][] Xsi_s,
                                      double[] extT, double[] t0, double[][] P0, double dx, int maxIt, double tol) {

        // get number of types and time steps
        int n = cdf[0].length;
        int m = t0.length;

        // initialize matrix
        double[][] X0 = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int x = 0; x < n; x++) {
                X0[i][x] = rho * (1 - cdf[i][x]);
            }
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[][] X = X0;

        // iterate
        while (err > tol && it < maxIt) {
            double[][] Xi = new double[m][n];

            for (int j = 0; j < n; j++) {
                // get vectors for convolution
                double[] y = new double[m];
                for (int i = 0; i < m; i++) { // multiply elementwise on times
                    for (int k = 0; k < n; k++) { // sum over all types k
                        y[i] += 2 * Xsi_as[j][k] * (X[i][j] * P0[i][k] + X[i][k] * P0[i][j])
                                + Xsi_s[j][k] * X[i][k] * P0[i][k];
                    }
                }
                // extract column from the pdf matrix
                Complex[] Ft = new Complex[m*2];
                for (int i = 0; i < m*2; i++) {
                    Ft[i] = pdfFFT[i][j];
                }

                // partially convolve
                double[] I = convolveFFT(Ft, y, m, dx);

                // sum
                for (int i = 0; i < m; i++) {
                    Xi[i][j] = X0[i][j] + 2 * (1 - d[j]) * I[i];
                }
            }

            // compute error
            err = getMatrixError(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxIt) {
            System.err.printf("calcP1 Warning: max iterations reached with error: %.2f%n", err);
        }

        // interpolate: evaluate at tip times
        int ntips = extT.length;
        double[][] P1 = new double[ntips][n];

        for (int j = 0; j < n; j++) {
            // extract column from X
            double[] Xj = new double[m];
            for (int i = 0; i < m; i++) {
                Xj[i] = X[i][j];
            }
            UnivariateFunction function = interpolator.interpolate(t0, Xj);

            for (int i = 0; i < ntips; i++) {
                P1[i][j] = function.value(extT[i]);
            }
        }

        return P1;
    }


    // Function for calculating branch probabilities
    public static double[][] calcMTB(double[] a, double[] b, double[] d, double[][] Xsi_as, double[][] Xsi_s,
                                     double[] s, double[] e, double[] t0, double[][] P0,
                                     int maxIt, double tol, int m) {

        // get number of types and branches
        int ntypes = a.length;
        int nbranches = s.length;
        assert e.length == nbranches; // for each branch, a start and end time must be given

        // interpolate P0
        List<UnivariateFunction> functions = new ArrayList<>();
        IntStream.range(0, ntypes)
                .parallel() // for each type j
                .forEach(j -> {
                    // extract column from P0
                    double[] P0j = new double[t0.length];
                    for (int i = 0; i < t0.length; i++) {
                        P0j[i] = P0[i][j];
                    }
                    UnivariateFunction function = interpolator.interpolate(t0, P0j);
                    functions.add(function);
                });

        double[][] B = new double[nbranches][ntypes];

        // for each branch:
        IntStream.range(0, nbranches)
                .parallel()
                .forEach(x -> {
                    double sx = s[x]; // start of the branch
                    double ex = e[x]; // end of the branch

                    // generate linearly spaced values between start and end
                    double[] tSeq = new double[m];
                    double[] age_seq = new double[m];
                    double dx = (ex - sx) / m;
                    for (int i = 0; i < m; i++) {
                        tSeq[i] = sx + dx * (i + 1);
                        age_seq[i] = tSeq[i] - sx;
                    }
                    assert tSeq[m - 1] == ex;

                    // calculate the FFT from PDF of the gamma distribution and interpolate P0 per type, initialize matrix
                    Complex[][] pdfFFT = new Complex[m*2][ntypes];
                    double[][] P = new double[m][ntypes];
                    double[][] X0 = new double[m][ntypes];

                    for (int j = 0; j < ntypes; j++) {
                        double[] pdf = new double[m];
                        GammaDistribution gammaDist = new GammaDistribution(b[j], a[j]);
                        for (int i = 0; i < m; i++) {
                            pdf[i] = gammaDist.density(age_seq[i]); // get density
                            P[i][j] = functions.get(j).value(tSeq[i]); // interpolate P0
                            X0[i][j] = (1 - d[j]) * pdf[i]; // initalize matrix
                        }

                        // perform FFT
                        Complex[] Ft = fft.transform(padZeros(pdf), TransformType.FORWARD);
                        for (int i = 0; i < m*2; i++) {
                            pdfFFT[i][j] = Ft[i];
                        }
                    }

                    // set up iteration
                    double err = 1;
                    int it = 0;
                    double[][] X = X0;

                    // iterate
                    while (err > tol && it < maxIt) {
                        double[][] Xi = new double[m][ntypes];

                        for (int j = 0; j < ntypes; j++) {

                            // get vectors for convolution
                            double[] y = new double[m];
                            for (int k = 0; k < ntypes; k++) { // sum over all types k
                                for (int i = 0; i < m; i++) { // multiply elementwise on times
                                    y[i] += Xsi_as[j][k] * (X[i][j] * P[i][k] + P[i][j] * X[i][k]) +
                                            Xsi_s[j][k] * X[i][k] * P[i][k];
                                }
                            }

                            // extract column from the pdf matrix
                            Complex[] Ft = new Complex[m * 2];
                            for (int i = 0; i < m * 2; i++) {
                                Ft[i] = pdfFFT[i][j];
                            }

                            // partially convolve
                            double[] I = convolveFFT(Ft, y, m, dx);

                            // sum
                            for (int i = 0; i < m; i++) {
                                Xi[i][j] = X0[i][j] + 2 * (1 - d[j]) * I[i];
                            }
                        }

                        // compute error
                        err = getMatrixError(X, Xi);

                        // update
                        X = Xi;
                        it++;
                    }

                    if (it == maxIt) {
                        System.err.printf("calcB branch %d Warning: max iterations reached with error: %.2f%n", x, err);
                    }

                    // take last row only
                    B[x] = X[m - 1];
                });

        return B;
    }


    // Helper method for calculating the 1-(Manhattan)-distance between two matrices,
    // which is the maximum absolute column sum of the matrices substracted element-wise
    private static double getMatrixError(double[][] X, double[][] Y) {
        int m = X.length; // number of rows
        int n = X[0].length;  // number of columns

        double maxColSum = 0;

        // iterate through each column
        for (int col = 0; col < n; col++) {
            double colSum = 0;

            // sum the absolute differences for this column
            for (int row = 0; row < m; row++) {
                colSum += Math.abs(X[row][col] - Y[row][col]);
            }

            // Update maxColumnSum if this column has a larger sum
            if (colSum > maxColSum) {
               maxColSum = colSum;
            }
        }

        return maxColSum;
    }

}




