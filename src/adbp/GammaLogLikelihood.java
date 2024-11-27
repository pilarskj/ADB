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

    private static FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
    private static EuclideanDistance norm = new EuclideanDistance();
    private static LinearInterpolator interpolator = new LinearInterpolator();


    // Main function for calculating the likelihood of the tree
    /* Inputs
    a: scale, b: shape, d: death probability, rho: sampling probability, origin: age of the tree
    intS: start times of internal branches, intE: end times of internal branches, extE: end times of external branches (backwards in time)
    Options for solving integral equations -
    maxIt: maximum number of iterations, tolP: error tolerance for P0 and P1, tolB: error tolerance for branch probabilities (B),
    mP: step size for P0 and P1, mB: step size for B, approx: approximate B?
     */
    public static double calcLogLikelihood(double a, int b, double d, double rho, double origin,
                                           double[] intS, double[] intE, double[] extE,
                                           int maxIt, double tolP, double tolB, int mP, int mB, boolean approx) {
        // initialize distribution
        GammaDistribution gammaDist = new GammaDistribution(b, a);

        // generate linearly spaced values between 0 and origin
        int m = mP; // the number of time steps must be a power of 2 (required by FFT!)

        double[] tSeq = new double[m];
        final double dx = origin / m;
        for (int i = 0; i < m; i++) {
            tSeq[i] = dx * (i + 1);
        }
        assert tSeq[m - 1] == origin;

        // calculate CDF and FFT from PDF from 0 to origin
        double[] pdf = new double[m];
        double[] cdf = new double[m];
        for (int i = 0; i < m; i++) {
            pdf[i] = gammaDist.density(tSeq[i]);
            cdf[i] = gammaDist.cumulativeProbability(tSeq[i]);
        }
        Complex[] pdfFFT = fft.transform(padZeros(pdf), TransformType.FORWARD);

        // calculate extinction probability over time
        final double[] P0 = calcP0(pdfFFT, cdf, d, rho, dx, maxIt, tolP);

        // calculate probability of single descendants at tips
        final double[] P1 = calcP1(pdfFFT, cdf, P0, d, rho, extE, tSeq, dx, maxIt, tolP);

        // calculate probabilities of internal branches
        double[] B;
        if (approx) {
            B = approxB(a, b, d, intS, intE, tSeq, P0, tolB);
        } else {
            B = calcB(a, b, d, intS, intE, tSeq, P0, maxIt, tolB, mB);
        }

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
    public static double[] calcP0(Complex[] pdfFFT, double[] cdf, double d, double rho, double dx, int maxIt, double tol) {

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
            double[] I = convolveFFT(pdfFFT, y, n, dx);

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
    public static double[] calcP1(Complex[] pdfFFT, double[] cdf, double[] P0, double d, double rho,
                                  double[] extT, double[] tSeq, double dx, int maxIt, double tol) {

        // get length
        int n = tSeq.length;

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
            double[] I = convolveFFT(pdfFFT, y, n, dx);

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
        UnivariateFunction function = interpolator.interpolate(tSeq, X);
        double[] P1 = new double[extT.length];
        for (int i = 0; i < extT.length; i++) {
            P1[i] = function.value(extT[i]);
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
                    double[] tSeq = new double[m];
                    double[] age_seq = new double[m];
                    double dx = (ex - sx) / m;
                    for (int i = 0; i < m; i++) {
                        tSeq[i] = sx + dx * (i + 1);
                        age_seq[i] = tSeq[i] - sx;
                    }
                    assert age_seq[m - 1] == ex - sx;

                    // calculate the PDF of the gamma distribution and compute P0
                    double[] pdf = new double[m];
                    double[] P = new double[m];
                    for (int i = 0; i < m; i++) {
                        pdf[i] = gammaDist.density(age_seq[i]);
                        P[i] = function.value(tSeq[i]);
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
        // Ugne: you could use cashing to avoid creating a new GammaDistribution object for each iteration.
        // Basically create a map where you store the distributions for each i and reuse them.
        // If you know the maximum i you will need, you can precompute them all and store them in a list.
        // If you don't know the maximum i, whenever you need gammaDist for particular i, you can use a map and check
        // if the distribution for i is already computed. And if not, you compute it and store it in the map.
        
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

                    // Ugne: Below you could use binary search instead of loop in the findClosestIndex method.
                    // This would make the complexity of the method log(n) instead of n.
                    // But you need to sort t0 array first, which would be an additional step.
                    // If you have a lot of branches and the t0 array is long, it may be worth it.
                    // Then you can sort it once, before repeatedly using in this stream.
                    // Be careful to check if sorting this array does not affect some dependent calculations.
                    // Otherwise, if looping is actually faster, you may change the method to find the closest index for two doubles at the same time.
                    // Since you anyway need to loop through all t0, you could find both indices for sx and ex at the same time.

                    // get average P0 over branch
                    double[] P0Slice = Arrays.copyOfRange(P0, findClosestIndex(t0, sx), findClosestIndex(t0, ex) + 1);
                    double P0M = getMean(P0Slice);

                    // initialize the approximation
                    int k = (int) ((ex - sx) / (a * b)); // for this k, the term b_k will be maximal
                    double branch_prob = 0;

                    // int num_terms = 0;
                    int i = k; // increase k
                    double term = 1;
                    while (term > tol) {
                        GammaDistribution gammaDist = new GammaDistribution((i+1) * b, a);
                        term = Math.pow(2, i) * Math.pow(1-d, i+1) * Math.pow(P0M, i) * gammaDist.density(ex - sx); // calculate b_i
                        branch_prob += term;
                        i++;
                        // num_terms++;
                    }
                    int j = k-1; // decrease k
                    term = 1;
                    while (term > tol & j >= 0) {
                        GammaDistribution gammaDist = new GammaDistribution((j+1) * b, a);
                        term = Math.pow(2, j) * Math.pow(1-d, j+1) * Math.pow(P0M, j) * gammaDist.density(ex - sx); // calculate b_j
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




