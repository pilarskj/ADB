package adbp;

import beast.base.core.Log;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;
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
    public static double calcLogLikelihood(double a, Number b, double d, double rho, double origin,
                                           double[] intS, double[] intE, double[] extE,
                                           int maxIt, double tolP, double tolB, int mP, int mB, boolean approx) {

        /*
        // use much simpler calculation for BD case
        if (b.doubleValue() == 1 && d != 0.5) {
            return calcBDLogLikelihood(a, d, rho, origin, intS, extE);
        } */


        // initialize distribution
        GammaDistribution gammaDist = new GammaDistribution(b.doubleValue(), a);

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
            pdf[i] = Math.exp(gammaDist.logDensity(tSeq[i])); // gammaDist.density(tSeq[i]) - use log to prevent underflow
            cdf[i] = gammaDist.cumulativeProbability(tSeq[i]);
        }
        Complex[] pdfFFT = fft.transform(padZeros(pdf), TransformType.FORWARD);

        // calculate extinction probability over time
        final double[] P0 = calcP0(pdfFFT, cdf, d, rho, dx, maxIt, tolP);
        /* // print error if P0 is not in range
        if (Arrays.stream(P0).min().getAsDouble() < 0 || Arrays.stream(P0).max().getAsDouble() > 1) {
            Log.debug.println("P0 not in range [0,1]");
        } */

        // calculate probability of single descendants at tips
        final double[] P1 = calcP1(pdfFFT, cdf, P0, d, rho, extE, tSeq, dx, maxIt, tolP);
        /* // print error if P1 is not in range
        if (Arrays.stream(P1).min().getAsDouble() < 0 || Arrays.stream(P1).max().getAsDouble() > 1) {
            Log.debug.println("P1 not in range [0,1]");
        } */

        // calculate probabilities of internal branches
        double[] B;
        if (approx) {
            B = approxB(a, b.intValue(), d, intS, intE, tSeq, P0, tolB);
        } else {
            // extend arrays (to calculate probabilities of very short internal branches)
            double[] extSeq = new double[tSeq.length + 1];
            double[] extP0 = new double[P0.length + 1];
            extSeq[0] = 0;
            extP0[0] = 1 - rho;
            System.arraycopy(tSeq, 0, extSeq, 1, tSeq.length);
            System.arraycopy(P0, 0, extP0, 1, P0.length);
            B = calcB(a, b, d, intS, intE, extSeq, extP0, maxIt, tolB, mB);
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

        // sum (with conditioning on survival)
        return -Math.log(1 - P0[m - 1]) + logP1 + logB;
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

        // for asserting convergence
        if (it == maxIt) {
            // print error in debug mode
            Log.debug.println("calcP0: max iterations reached with error:  " + err);
        }

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

        if (it == maxIt) {
            Log.debug.println("calcP1: max iterations reached with error:  " + err);
        }

        // interpolate
        // extend arrays (to calculate probabilities of very short external branches)
        double[] extSeq = new double[tSeq.length + 1];
        double[] extX = new double[X.length + 1];
        extSeq[0] = 0;
        extX[0] = rho;
        System.arraycopy(tSeq, 0, extSeq, 1, tSeq.length);
        System.arraycopy(X, 0, extX, 1, X.length);

        UnivariateFunction function = interpolator.interpolate(extSeq, extX);
        double[] P1 = new double[extT.length];
        for (int i = 0; i < extT.length; i++) {
            P1[i] = function.value(extT[i]); // interpolate
        }

        return P1;
    }


    // Function for calculating branch probabilities
    public static double[] calcB(double a, Number b, double d,
                                 double[] s, double[] e, double[] t0, double[] P0,
                                 int maxIt, double tol, int m) {

        // initialize distribution
        GammaDistribution gammaDist = new GammaDistribution(b.doubleValue(), a);

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
                        pdf[i] = Math.exp(gammaDist.logDensity(age_seq[i]));
                        P[i] = function.value(tSeq[i]);
                    }
                    double maxP = Arrays.stream(P).max().getAsDouble();
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
                        err = norm.compute(X, Xi);

                        // update
                        X = Xi;
                        it++;
                    }

                    if (it == maxIt) {
                        Log.debug.println("calcB: max iterations reached with error:  " + err);
                    }

                    // take last element
                    B[x] = X[m - 1];
                });

        return B;
    }


    // Function for approximating branch probabilities
    public static double[] approxB(double a, int b, double d,
                                   double[] s, double[] e, double[] t0, double[] P0,
                                   double tol) {

        // use cashing to avoid creating a new GammaDistribution object with shape i*b for i=1,2,... in each iteration
        // create a map to store the distributions for each i (once needed) and re-use them
        // use a thread-safe and dynamic map
        ConcurrentHashMap<Integer, GammaDistribution> gammaCache = new ConcurrentHashMap<>();

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
                    // use binary search to find closest indices to the branch lengths in t0 (t0 is sorted per definition!)
                    double[] P0Slice = Arrays.copyOfRange(P0, findClosestIndex(t0, sx), findClosestIndex(t0, ex) + 1);
                    double P0M = getMean(P0Slice);

                    // initialize the approximation
                    int k = (int) ((ex - sx) / (a * b)); // for this k, the term b_k will be maximal
                    double branchProb = 0;

                    // int num_terms = 0;
                    int i = k; // increase k
                    double term = 1;
                    while (term > tol) {
                        GammaDistribution gammaDist;
                        // check if GammaDistribution with the given shape is already in the cache
                        if (gammaCache.containsKey(i + 1)) {
                            gammaDist = gammaCache.get(i + 1);
                        } else { // otherwise, add it
                            gammaDist = new GammaDistribution((i + 1) * b, a);
                            gammaCache.put(i + 1, gammaDist);
                        }
                        term = Math.pow(2, i) * Math.pow(1 - d, i + 1) * Math.pow(P0M, i) * Math.exp(gammaDist.logDensity(ex - sx)); // calculate b_i
                        branchProb += term;
                        i++;
                    }
                    int j = k - 1; // decrease k
                    term = 1;
                    while (term > tol & j >= 0) {
                        GammaDistribution gammaDist;
                        // check if GammaDistribution with the given shape is already in the cache
                        if (gammaCache.containsKey(j + 1)) {
                            gammaDist = gammaCache.get(j + 1);
                        } else { // otherwise, add it
                            gammaDist = new GammaDistribution((j + 1) * b, a);
                            gammaCache.put(j + 1, gammaDist);
                        }
                        term = Math.pow(2, j) * Math.pow(1 - d, j + 1) * Math.pow(P0M, j) * Math.exp(gammaDist.logDensity(ex - sx)); // calculate b_j
                        branchProb += term;
                        j--;
                    }

                    B[x] = branchProb;
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


    // Helper method to find the index of the closest value in a sorted array using binary search
    // https://stackoverflow.com/questions/30245166/find-the-nearest-closest-value-in-a-sorted-list
    // complexity log(n) instead of n in a loop
    public static int findClosestIndex(double[] array, double value) {
        // if value is at boundaries
        if (value <= array[0]) {
            return 0;
        }
        if (value >= array[array.length - 1]) {
            return array.length - 1;
        }

        // do binary search
        int index = Arrays.binarySearch(array, value);

        if (index >= 0) { // exact match found
            return index;
        } else { // no exact match: binarySearch returns (-(insertion point) - 1)
            int insertionPoint = -(index + 1);

            // return the index of the closest value
            if ((value - array[insertionPoint - 1]) <= (array[insertionPoint] - value)) { // value is closer to the previous value
                return insertionPoint - 1;
            } else { // value is closer to the next value
                return insertionPoint;
            }
        }
    }


    // Helper method to calculate the mean of an array
    private static double getMean(double[] array) {
        double sum = 0;
        for (double num : array) {
            sum += num;
        }
        return sum / array.length;
    }


    // simpler function for calculating the log likelihood of a birth-death-sampling tree (Stadler, JTB 2010) - for testing
    public static double calcBDLogLikelihood(double a, double d, double rho, double origin,
                                             double[] intS, double[] extE) {

        // get birth and death rate
        double lambda = (1 - d) / a;
        double mu = d / a;

        // get all bifurcation times
        double[] t = new double[intS.length + 1];
        System.arraycopy(intS, 0, t, 0, intS.length);
        t[t.length - 1] = origin;

        // get branch probabilities (no sampling through time)
        double c1 = lambda - mu;
        double c2 = -(lambda - mu - 2 * lambda * rho) / c1;
        double[] B = new double[t.length];
        IntStream.range(0, t.length)
                .parallel()
                .forEach(i -> {
                    B[i] = 2 * (1 - Math.pow(c2, 2)) +
                            Math.exp(-c1 * t[i]) * Math.pow(1 - c2, 2) +
                            Math.exp(c1 * t[i]) * Math.pow(1 + c2, 2);
                });

        // make log and sum
        double logB = 0;
        for (int i = 0; i < B.length; i++) {
            logB += Math.log(B[i]);
        }

        // condition on survival
        double p0 = 1 - (rho * (lambda - mu)) / (rho * lambda + (lambda * (1 - rho) - mu) * Math.exp(-(lambda - mu) * origin));

        // get log likelihood
        int n = extE.length; // number of tips
        double logL = -Math.log(1 - p0) + (n - 1) * Math.log(lambda) + n * Math.log(4 * rho) - logB;

        return logL;
    }
}