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


public class GammaLogLikelihood {

    public static double calcLogLikelihood(double rho, double a, double b, double d,
                                           double t_or, double[] int_s, double[] int_e, double[] ext_e,
                                           int mP, int mB, int maxit) {
        // set large m for P0 and P1
        int m = mP;
        // m must be a power of 2 (required by FFT!)
        // generate linearly spaced values between 0 and origin
        double[] t_seq = new double[m];
        double dx = t_or / m; // (m - 1)
        for (int i = 0; i < m; i++) {
            t_seq[i] = dx * (i + 1); // i
        }
        assert t_seq[m - 1] == t_or;

        // calculate extinction probability over time
        double[] P0 = calcP0(rho, a, b, d, t_seq, dx, maxit);

        // calculate probability of single descendants at tips
        double[] P1 = calcP1(rho, a, b, d, ext_e, t_seq, P0, dx, maxit);

        // reduce m for B
        m = mB;

        // calculate probabilities of internal branches
        // double[] B = calcB(a, b, d, int_s, int_e, t_seq, P0, m, maxit);
        double[] B = approxB(a, b, d, int_s, int_e, t_seq, P0);

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
        double logL = logP1 + logB;
        return logL;
    }


    // Function for calculating the extinction probability
    public static double[] calcP0(double rho, double a, double b, double d, double[] t, double dx, int maxit) {

        // get length
        int n = t.length;

        // calculate the PDF and CDF of the gamma distribution
        double[] pdf = new double[n];
        double[] cdf = new double[n];
        GammaDistribution gammaDist = new GammaDistribution(b, a);
        for (int i = 0; i < n; i++) {
            pdf[i] = gammaDist.density(t[i]);
            cdf[i] = gammaDist.cumulativeProbability(t[i]);
        }

        // initialize
        double[] X0 = new double[n];
        for (int i = 0; i < n; i++) {
            X0[i] = (1 - rho) * (1 - cdf[i]) + d * cdf[i];
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[] X = X0;

        // perform FFT
        double[] ft = padZeros(pdf);

        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] Ft = fft.transform(ft, TransformType.FORWARD);

        // iterate
        while (err > 1e-12 && it < maxit) {

            // square
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = X[i] * X[i];
            }

            // partially convolve
            double[] I = convolveFFT(Ft, y, n, dx);

            // sum
            double[] Xi = new double[n];
            for (int i = 0; i < n; i++) {
                Xi[i] = X0[i] + (1 - d) * I[i];
            }

            // compute error
            EuclideanDistance norm = new EuclideanDistance();
            err = norm.compute(Xi, X);

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
    public static double[] calcP1(double rho, double a, double b, double d,
                                  double[] t, double[] t0, double[] P0, double dx, int maxit) {

        // get length
        int n = t0.length;

        // calculate the PDF and CDF of the gamma distribution
        double[] pdf = new double[n];
        double[] cdf = new double[n];
        GammaDistribution gammaDist = new GammaDistribution(b, a);
        for (int i = 0; i < n; i++) {
            pdf[i] = gammaDist.density(t0[i]);
            cdf[i] = gammaDist.cumulativeProbability(t0[i]);
        }

        // initialize
        double[] X0 = new double[n];
        for (int i = 0; i < n; i++) {
            X0[i] = rho * (1 - cdf[i]);
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[] X = X0;

        // perform FFT
        double[] ft = padZeros(pdf);
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex [] Ft = fft.transform(ft, TransformType.FORWARD);

        // iterate
        while (err > 1e-12 && it < maxit) {

            // multiply
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = P0[i] * X[i];
            }

            // partially convolve
            double[] I = convolveFFT(Ft, y, n, dx);

            // sum
            double[] Xi = new double[n];
            for (int i = 0; i < n; i++) {
                Xi[i] = X0[i] + 2 * (1 - d) * I[i];
            }

            // compute error
            EuclideanDistance norm = new EuclideanDistance();
            err = norm.compute(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxit) {
            System.err.printf("calcP1 Warning: max iterations reached with error: %.2f%n", err);
        }

        // interpolate
        LinearInterpolator interpolator = new LinearInterpolator();
        UnivariateFunction function = interpolator.interpolate(t0, X);
        double[] P1 = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            P1[i] = function.value(t[i]);
        }
        return P1;
    }


    // Function for calculating branch probabilities
    public static double[] calcB(double a, double b, double d,
                                 double[] s, double[] e, double[] t0, double[] P0,
                                 int m, int maxit) {

        // interpolate P0
        LinearInterpolator interpolator = new LinearInterpolator();
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
                    double dx = (ex - sx) / m; // (m - 1)
                    for (int i = 0; i < m; i++) {
                        t_seq[i] = sx + dx * (i + 1); // i
                        age_seq[i] = t_seq[i] - sx;
                    }
                    assert t_seq[m - 1] == ex;

                    // calculate the PDF of the gamma distribution
                    double[] pdf = new double[m];
                    GammaDistribution gammaDist = new GammaDistribution(b, a);
                    for (int i = 0; i < m; i++) {
                        pdf[i] = gammaDist.density(age_seq[i]);
                    }

                    // interpolate
                    double[] P = new double[m];
                    for (int i = 0; i < m; i++) {
                        P[i] = function.value(t_seq[i]);
                    }

                    // initialize
                    double[] X0 = new double[m];
                    for (int i = 0; i < m; i++) {
                        X0[i] = (1 - d) * pdf[i];
                    }

                    // set up iteration
                    double err = 1;
                    int it = 0;
                    double[] X = X0;

                    // perform FFT
                    double[] ft = padZeros(pdf);
                    FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
                    Complex [] Ft = fft.transform(ft, TransformType.FORWARD);

                    // iterate
                    while (err > 1e-12 && it < maxit) {

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

                    if (it == maxit) {
                        System.err.printf("calcB branch %d Warning: max iterations reached with error: %.2f%n", x, err);
                    }

                    // take last element
                    B[x] = X[m - 1];
                });

        return B;
    }



    // Function for approximating branch probabilities
    public static double[] approxB(double a, double b, double d,
                                 double[] s, double[] e, double[] t0, double[] P0) {

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
                    GammaDistribution gammaDist = new GammaDistribution(b, a);
                    double branch_prob = (1 - d) * gammaDist.density(ex - sx);
                    double min_num_terms = (ex - sx) / (a * b);

                    int i = 1;
                    double update = 1;
                    while (update > 1e-6 || i <= min_num_terms) {
                        GammaDistribution newGammaDist = new GammaDistribution((i+1) * b, a);
                        update = Math.pow(2, i) * Math.pow(1-d, i+1) * Math.pow(P0_bar, i) * newGammaDist.density(ex - sx);
                        branch_prob += update;
                        i++;
                    }

                    // take last element
                    B[x] = branch_prob;
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


    /* redundant calculation for P0 and P1
    private static Complex[] getPDF(double[] t, double a, double b) {
        int n = t.length;

        // calculate the PDF of the gamma distribution
        double[] pdf = new double[n];
        GammaDistribution gammaDist = new GammaDistribution(b, a);
        for (int i = 0; i < n; i++) {
            pdf[i] = gammaDist.density(t[i]);
        }

        // perform FFT
        double[] ft = padZeros(pdf);
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex [] Ft = fft.transform(ft, TransformType.FORWARD);

        return Ft;
    }
    */
}




