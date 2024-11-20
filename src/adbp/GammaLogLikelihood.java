package adbp;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.Arrays;
import java.util.stream.IntStream;

/*
Class for solving the equations for P0, P1 and B
and for calculating the likelihood of a tree based on the parameters and branching times
 */
public class GammaLogLikelihood {

    // can those be initialized once and re-used?
    public static FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
    public static EuclideanDistance norm = new EuclideanDistance();
    public static LinearInterpolator interpolator = new LinearInterpolator();

    // Main function for calculating the likelihood of the tree
    public static double calcLogLikelihood(double a, int b, double d, double rho, double origin,
                                           double[] int_s, double[] int_e, double[] ext_e,
                                           int maxIt, double tolP, double tolB, int mP, int mB) {
        // initialize distribution
        GammaDistribution gammaDist = new GammaDistribution(b, a);

        // generate linearly spaced values between 0 and origin
        int m = mP; // the number of time steps must be a power of 2 (required by FFT!)
        double[] t_seq = new double[m];
        double dx = origin / m;
        for (int i = 0; i < m; i++) {
            t_seq[i] = dx * (i + 1);
        }
        assert t_seq[m - 1] == origin;

        // calculate CDF and FFT from PDF from 0 to origin
        double[] pdf = new double[m];
        double[] cdf = new double[m];
        for (int i = 0; i < m; i++) {
            pdf[i] = gammaDist.density(t_seq[i]);
            cdf[i] = gammaDist.cumulativeProbability(t_seq[i]);
        }
        Complex[] pdf_FFT = fft.transform(padZeros(pdf), TransformType.FORWARD);

        // calculate extinction probability over time
        double[] P0 = calcP0(pdf_FFT, cdf, d, rho, dx, maxIt, tolP);

        // calculate probability of single descendants at tips
        double[] P1 = calcP1(pdf_FFT, cdf, P0, d, rho, ext_e, t_seq, dx, maxIt, tolP);

        // calculate probabilities of internal branches
        // double[] B = calcB(a, b, d, int_s, int_e, t_seq, P0, maxIt, tolB, mB);
        double[] B = approxB(a, b, d, int_s, int_e, t_seq, P0, tolB);

        // make log and sum
        double logP1 = 0;
        for (int i = 0; i < P1.length; i++) {
            logP1 += Math.log(P1[i]);
        }
        double logB = 0;
        for (int i = 0; i < B.length; i++) {
            logB += Math.log(B[i]);
        }

        // sum
        return logP1 + logB;
    }


    // Function for calculating the extinction probability
    public static double[] calcP0(Complex[] pdf_FFT, double[] cdf, double d, double rho, double dx, int maxIt, double tol) {

        // get length
        int n = cdf.length;

        // initialize
        double[] X0 = new double[n];
        for (int i = 0; i < n; i++) {
            X0[i] = (1 - rho) * (1 - cdf[i]) + d * cdf[i];
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[] X = X0;

        // iterate
        while (err > tol && it < maxIt) {

            // square
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = X[i] * X[i];
            }

            // partially convolve
            double[] I = convolveFFT(pdf_FFT, y, n, dx);

            // sum
            double[] Xi = new double[n];
            for (int i = 0; i < n; i++) {
                Xi[i] = X0[i] + (1 - d) * I[i];
            }

            // compute error
            err = norm.compute(Xi, X);

            // update
            X = Xi;
            it++;
        }

        /* for asserting convergence
        if (it == maxIt) {
            System.err.printf("calcP0 Warning: max iterations reached with error: %.2f%n", err); // check if error is huge!
        }
         */

        return X;
    }


    // Function for calculating the probability of a single descendant
    public static double[] calcP1(Complex[] pdf_FFT, double[] cdf, double[] P0, double d, double rho,
                                  double[] ext_t, double[] t_seq, double dx, int maxIt, double tol) {

        // get length
        int n = t_seq.length;

        // initialize
        double[] X0 = new double[n];
        for (int i = 0; i < n; i++) {
            X0[i] = rho * (1 - cdf[i]);
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[] X = X0;

        // iterate
        while (err > tol && it < maxIt) {

            // multiply
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = P0[i] * X[i];
            }

            // partially convolve
            double[] I = convolveFFT(pdf_FFT, y, n, dx);

            // sum
            double[] Xi = new double[n];
            for (int i = 0; i < n; i++) {
                Xi[i] = X0[i] + 2 * (1 - d) * I[i];
            }

            // compute error
            err = norm.compute(X, Xi);

            // update
            X = Xi;
            it++;
        }

        /*
        if (it == maxIt) {
            System.err.printf("calcP1 Warning: max iterations reached with error: %.2f%n", err);
        }
         */

        // interpolate
        UnivariateFunction function = interpolator.interpolate(t_seq, X);
        double[] P1 = new double[ext_t.length];
        for (int i = 0; i < ext_t.length; i++) {
            P1[i] = function.value(ext_t[i]);
        }

        return P1;
    }


    // Function for calculating branch probabilities
    public static double[] calcB(double a, int b, double d,
                                 double[] s, double[] e, double[] t0, double[] P0,
                                 int maxIt, double tol, int m) {

        // initialize distribution
        GammaDistribution gammaDist = new GammaDistribution(b, a);

        // interpolate P0
        UnivariateFunction function = interpolator.interpolate(t0, P0);

        // get number of branches
        assert s.length == e.length; // for each branch, a start and end time must be given
        int n = s.length;

        double[] B = new double[n];
        // for each branch:
        // use Stream for parallelization
        IntStream.range(0, n)
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
                    assert age_seq[m - 1] == ex - sx;

                    // calculate the PDF of the gamma distribution and compute P0
                    double[] pdf = new double[m];
                    double[] P = new double[m];
                    for (int i = 0; i < m; i++) {
                        pdf[i] = gammaDist.density(age_seq[i]);
                        P[i] = function.value(t_seq[i]);
                    }
                    Complex[] Ft = fft.transform(padZeros(pdf), TransformType.FORWARD); // perform FFT

                    // initialize
                    double[] X0 = new double[m];
                    for (int i = 0; i < m; i++) {
                        X0[i] = (1 - d) * pdf[i];
                    }

                    // set up iteration
                    double err = 1;
                    int it = 0;
                    double[] X = X0;

                    // iterate
                    while (err > tol && it < maxIt) {

                        // multiply
                        double[] y = new double[m];
                        for (int i = 0; i < m; i++) {
                            y[i] = P[i] * X[i];
                        }

                        // partially convolve
                        double[] I = convolveFFT(Ft, y, m, dx);

                        // sum
                        double[] Xi = new double[m];
                        for (int i = 0; i < m; i++) {
                            Xi[i] = X0[i] + 2 * (1 - d) * I[i];
                        }

                        // compute error
                        EuclideanDistance norm = new EuclideanDistance();
                        err = norm.compute(X, Xi);

                        // update
                        X = Xi;
                        it++;
                    }

                    /*
                    if (it == maxIt) {
                        System.err.printf("calcB branch %d Warning: max iterations reached with error: %.2f%n", x, err);
                    }
                     */

                    // take last element
                    B[x] = X[m - 1];
                });

        return B;
    }


    // Function for approximating branch probabilities
    public static double[] approxB(double a, int b, double d,
                                   double[] s, double[] e, double[] t0, double[] P0,
                                   double tol) {

        // to check: initialize GammaDistribution with shape i*b for i=1,2,... here?

        // get number of branches
        assert s.length == e.length; // for each branch, a start and end time must be given
        int n = s.length;

        double[] B = new double[n];
        // for each branch:
        // use Stream for parallelization
        IntStream.range(0, n)
                .parallel()
                .forEach(x -> {
                    double sx = s[x]; // start of the branch
                    double ex = e[x]; // end of the branch

                    // get average P0 over branch
                    double[] P0_slice = Arrays.copyOfRange(P0, findClosestIndex(t0, sx), findClosestIndex(t0, ex) + 1);
                    double P0_bar = getMean(P0_slice);

                    // initialize the approximation
                    int k = (int) ((ex - sx) / (a * b)); // for this k, the term b_k will be maximal
                    double branch_prob = 0;

                    // int num_terms = 0;
                    int i = k; // increase k
                    double term = 1;
                    while (term > tol) {
                        GammaDistribution gammaDist = new GammaDistribution((i+1) * b, a);
                        term = Math.pow(2, i) * Math.pow(1-d, i+1) * Math.pow(P0_bar, i) * gammaDist.density(ex - sx); // calculate b_i
                        branch_prob += term;
                        i++;
                        // num_terms++;
                    }
                    int j = k-1; // decrease k
                    term = 1;
                    while (term > tol & j >= 0) {
                        GammaDistribution gammaDist = new GammaDistribution((j+1) * b, a);
                        term = Math.pow(2, j) * Math.pow(1-d, j+1) * Math.pow(P0_bar, j) * gammaDist.density(ex - sx); // calculate b_j
                        branch_prob += term;
                        j--;
                        // num_terms++;
                    }
                    // System.out.println("branch start " + sx + " end " + ex + " terms " + num_terms);
                    B[x] = branch_prob;
                });

        return B;
    }


    // Function for partial convolution using FFT
    public static double[] convolveFFT(Complex[] fx, double[] y, int n, double eps) {

        // perform FFT on padded y
        Complex[] fy = fft.transform(padZeros(y), TransformType.FORWARD);

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


    // Helper method to pad an array with 0 to its double length
    public static double[] padZeros(double[] x) {
        int n = x.length;
        double[] xp = new double[n * 2];
        System.arraycopy(x, 0, xp, 0, n);
        return xp;
    }


    // Helper method to find the index of the closest value in the array
    private static int findClosestIndex(double[] array, double value) {
        int closestIndex = 0;
        double minDiff = Math.abs(array[0] - value);
        for (int i = 1; i < array.length; i++) {
            double diff = Math.abs(array[i] - value);
            if (diff < minDiff) {
                minDiff = diff;
                closestIndex = i;
            }
        }
        return closestIndex;
    }


    // Helper method to calculate the mean of an array
    private static double getMean(double[] array) {
        double sum = 0;
        for (double num : array) {
            sum += num;
        }
        return sum / array.length;
    }
}




