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

public class MTLogLikelihood {

    public static double calcLogLikelihood(double[] a, double[] b, double[] d, double rho,
                                           double[][] Xsi_as, double[][] Xsi_s,
                                           double t_or, int type_or,
                                           double[] int_s, double[] int_e, double[] ext_e,
                                           double[] left_child, double[] right_child,
                                           int m, int mB, int maxit) {

        // get number of types
        int ntypes = a.length;
        assert b.length == ntypes;
        assert d.length == ntypes;

        // m must be a power of 2 (required by FFT!)
        // generate linearly spaced values between 0 and origin
        double[] t_seq = new double[m];
        double dx = t_or / m;
        for (int i = 0; i < m; i++) {
            t_seq[i] = dx * (i + 1);
        }
        assert t_seq[m - 1] == t_or;

        // calculate CDF and FFT from PDF from 0 to origin
        double[][] cdf = new double[m][ntypes];
        Complex[][] pdf_FFT = new Complex[m*2][ntypes];
        IntStream.range(0, ntypes)
                .parallel()
                .forEach(x -> {
                    // get density and cumulative probability
                    double[] pdf = new double[m];
                    GammaDistribution gammaDist = new GammaDistribution(b[x], a[x]);
                    for (int i = 0; i < m; i++) {
                        pdf[i] = gammaDist.density(t_seq[i]);
                        cdf[i][x] = gammaDist.cumulativeProbability(t_seq[i]);
                    }

                    // perform FFT
                    FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
                    Complex[] Ft = fft.transform(padZeros(pdf), TransformType.FORWARD);
                    for (int i = 0; i < m*2; i++) {
                        pdf_FFT[i][x] = Ft[i];
                    }
                });



        // calculate extinction probability over time
        double[][] P0 = calcP0(pdf_FFT, cdf, d, rho, Xsi_as, Xsi_s, t_seq, dx, maxit);

        // calculate probability of single descendants at tips
        double[][] P1 = calcP1(pdf_FFT, cdf, d, rho, Xsi_as, Xsi_s, ext_e, t_seq, P0, dx, maxit);

        // calculate probabilities of internal branches
        double[][] B = calcB(a, b, d,  Xsi_as, Xsi_s, int_s, int_e, t_seq, P0, mB, maxit);

        // start recursion
        double treeL = calcTreeLikelihood();

        return Math.log(treeL);
    }

    // Function for calculating the extinction probability
    public static double calcTreeLikelihood(double likelihood, int ntips, int ntypes,
                                            int stem_index, int stem_type, int[] left_child, int[] right_child,
                                            double[][] Xsi_as, double[][] Xsi_s,
                                            double[][] P1, double[][] B) {

        // at tips
        if (stem_index >= ntips - 1) {
            return P1[stem_index - (ntips - 1)][stem_type];
        } else {

            // recursion
            double subtreeL = 0;
            for (int j=0;  j < ntypes; j++){
                // get subtree (new stems)
                int right_stem = right_child[stem_index];
                int left_stem = left_child[stem_index];

                subtreeL += Xsi_as[stem_type][j] *
                        calcTreeLikelihood(likelihood, ntips, ntypes, right_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1, B) *
                        calcTreeLikelihood(likelihood, ntips, ntypes, left_stem, stem_type, left_child, right_child, Xsi_as, Xsi_s, P1, B) +
                        Xsi_as[stem_type][j] *
                                calcTreeLikelihood(likelihood, ntips, ntypes, right_stem, stem_type, left_child, right_child, Xsi_as, Xsi_s, P1, B) *
                                calcTreeLikelihood(likelihood, ntips, ntypes, left_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1, B) +
                        Xsi_s[stem_type][j] *
                                calcTreeLikelihood(likelihood, ntips, ntypes, right_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1, B) *
                                calcTreeLikelihood(likelihood, ntips, ntypes, left_stem, j, left_child, right_child, Xsi_as, Xsi_s, P1, B);
            }
            return B[stem_index][stem_type] * subtreeL;
        }
    }


    // Function for calculating the extinction probability
    public static double[][] calcP0(Complex[][] pdf_FFT, double[][] cdf, double[] d, double rho,
                                    double[][] Xsi_as, double[][] Xsi_s,
                                    double[] t, double dx, int maxit) {

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
        while (err > 1e-12 && it < maxit) {
            double[][] Xi = new double[m][n];

            IntStream.range(0, n)
                    .parallel() // for each type j
                    .forEach(j -> {
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
                            Ft[i] = pdf_FFT[i][j];
                        }

                        // partially convolve
                        double[] I = convolveFFT(Ft, y, m, dx);

                        // sum
                        for (int i = 0; i < m; i++) {
                            Xi[i][j] = X0[i][j] + (1 - d[j]) * I[i];
                        }
                    });

            // compute error
            err = getMatrixError(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxit) {
            System.err.printf("calcP0 Warning: max iterations reached with error: %.2f%n", err);
        }

        return X;
    }


    // Function for calculating the probability of a single descendant
    public static double[][] calcP1(Complex[][] pdf_FFT, double[][] cdf, double[] d, double rho,
                                    double[][] Xsi_as, double[][] Xsi_s,
                                    double[] t, double[] t0, double[][] P0, double dx, int maxit) {

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
        while (err > 1e-12 && it < maxit) {
            double[][] Xi = new double[m][n];

            IntStream.range(0, n)
                    .parallel() // for each type j
                    .forEach(j -> {
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
                            Ft[i] = pdf_FFT[i][j];
                        }

                        // partially convolve
                        double[] I = convolveFFT(Ft, y, m, dx);

                        // sum
                        for (int i = 0; i < m; i++) {
                            Xi[i][j] = X0[i][j] + 2 * (1 - d[j]) * I[i];
                        }
                    });

            // compute error
            err = getMatrixError(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxit) {
            System.err.printf("calcP1 Warning: max iterations reached with error: %.2f%n", err);
        }

        // interpolate: evaluate at tip times
        int ntips = t.length;
        double[][] P1 = new double[ntips][n];
        IntStream.range(0, n)
                .parallel() // for each type j
                .forEach(j -> {
                    // extract column from X
                    double[] Xj = new double[m];
                    for (int i = 0; i < m; i++) {
                        Xj[i] = X[i][j];
                    }
                    LinearInterpolator interpolator = new LinearInterpolator();
                    UnivariateFunction function = interpolator.interpolate(t0, Xj);

                    for (int i = 0; i < ntips; i++) {
                        P1[i][j] = function.value(t[i]);
                    }
                });

        return P1;
    }


    // Function for calculating branch probabilities
    public static double[][] calcB(double[] a, double[] b, double[] d, double[][] Xsi_as, double[][] Xsi_s,
                                   double[] s, double[] e, double[] t0, double[][] P0,
                                   int m, int maxit) {

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
                    LinearInterpolator interpolator = new LinearInterpolator();
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
                    double[] t_seq = new double[m];
                    double[] age_seq = new double[m];
                    double dx = (ex - sx) / m;
                    for (int i = 0; i < m; i++) {
                        t_seq[i] = sx + dx * (i + 1);
                        age_seq[i] = t_seq[i] - sx;
                    }
                    assert t_seq[m - 1] == ex;

                    // calculate the FFT from PDF of the gamma distribution and interpolate P0 per type, initialize matrix
                    Complex[][] pdf_FFT = new Complex[m*2][ntypes];
                    double[][] P = new double[m][ntypes];
                    double[][] X0 = new double[m][ntypes];

                    for (int j = 0; j < ntypes; j++) {
                        double[] pdf = new double[m];
                        GammaDistribution gammaDist = new GammaDistribution(b[x], a[x]);
                        for (int i = 0; i < m; i++) {
                            pdf[i] = gammaDist.density(age_seq[i]); // get density
                            P[i][j] = functions.get(j).value(t_seq[i]); // interpolate P0
                            X0[i][j] = (1 - d[j]) * pdf[i]; // initalize matrix
                        }

                        // perform FFT
                        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
                        Complex[] Ft = fft.transform(padZeros(pdf), TransformType.FORWARD);
                        for (int i = 0; i < m*2; i++) {
                            pdf_FFT[i][j] = Ft[i];
                        }
                    }

                    // set up iteration
                    double err = 1;
                    int it = 0;
                    double[][] X = X0;

                    // iterate
                    while (err > 1e-12 && it < maxit) {
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
                                Ft[i] = pdf_FFT[i][j];
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

                    if (it == maxit) {
                        System.err.printf("calcB branch %d Warning: max iterations reached with error: %.2f%n", x, err);
                    }

                    // take last row only
                    B[x] = X[m - 1];
                });

        return B;
    }


    // Function for partial convolution using FFT
    public static double[] convolveFFT(Complex[] fx, double[] y, int n, double eps) {

        // pad y
        double[] y_ext = padZeros(y);

        // perform FFT on y_ext
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] fy = fft.transform(y_ext, TransformType.FORWARD);

        // element-wise multiplication of fx and fy (convolution in Fourier space)
        Complex[] fz = new Complex[fx.length];
        for (int i = 0; i < fx.length; i++) {
            fz[i] = fx[i].multiply(fy[i]);
        }

        // perform inverse FFT to get the result back in time domain
        Complex[] z = fft.transform(fz, TransformType.INVERSE);

        // extract the real part and scale it by eps
        double[] z_real = new double[n];
        for (int i = 0; i < n; i++) {
            z_real[i] = z[i].getReal() * eps;
        }

        return z_real;
    }


    private static double[] padZeros(double[] x) {
        int n = x.length;

        // extend x
        double[] xp = new double[n * 2];
        System.arraycopy(x, 0, xp, 0, n);

        return xp;
    }


    // Function for calculating the 1-(Manhattan)-distance between two matrices,
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




