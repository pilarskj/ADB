import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.UnivariateFunction;

import static org.apache.commons.math3.complex.ComplexUtils.convertToComplex;

public class GammaLogLikelihood {

    public static double calcLogLikelihood(double rho, double a, double b,
                                           double t_or, double[] int_s, double[] int_e, double[] ext_e,
                                           int m, int maxit) {


        // generate linearly spaced values between 0 and origin
        double[] t_seq = new double[m];
        double dx = t_or / (m - 1);
        for (int i = 0; i < m; i++) {
            t_seq[i] = dx * i;
        }
        assert t_seq[m - 1] == t_or;

        // calculate extinction probability over time
        double[] P0 = calcP0(rho, a, b, t_seq, dx, maxit);

        // calculate probability of single descendants at tips
        double[] P1 = calcP1(rho, a, b, ext_e, t_seq, P0, dx, maxit);

        // calculate probabilities of internal branches
        double[] B = calcB(a, b, int_s, int_e, t_seq, P0, m, maxit);

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
    protected static double[] calcP0(double rho, double a, double b, double[] t, double dx, int maxit) {

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
            X0[i] = (1 - rho) * (1 - cdf[i]);
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[] X = X0;

        // perform FFT
        double[] ft = padReal(pdf);
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] Ft = fft.transform(convertToComplex(ft), TransformType.FORWARD);

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
                Xi[i] = X0[i] + I[i];
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
    protected static double[] calcP1(double rho, double a, double b, double[] t, double[] t0, double[] P0, double dx, int maxit) {

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
        double[] ft = padReal(pdf);
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex [] Ft = fft.transform(convertToComplex(ft), TransformType.FORWARD);

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
                Xi[i] = X0[i] + 2 * I[i];
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
    protected static double[] calcB(double a, double b,
                                    double[] s, double[] e, double[] t0, double[] P0,
                                    int m, int maxit) {

        // get number of branches
        assert s.length == e.length;
        int n = s.length;

        double[] B = new double[n];
        // for each branch:
        for (int x = 0; x < n; x++) {
            double sx = s[x]; // start of the branch
            double ex = e[x]; // end of the branch

            // generate linearly spaced values between start and end
            double[] t_seq = new double[m];
            double[] age_seq = new double[m];
            double dx = (ex - sx) / (m - 1);
            for (int i = 0; i < m; i++) {
                t_seq[i] = sx + dx * i;
                age_seq[i] = t_seq[i] - sx;
            }
            assert t_seq[m - 1] == ex;

            // calculate the PDF of the gamma distribution
            double[] pdf = new double[m];
            GammaDistribution gammaDist = new GammaDistribution(b, a);
            for (int i = 0; i < m; i++) {
                pdf[i] = gammaDist.density(age_seq[i]);
            }

            // interpolate P0
            LinearInterpolator interpolator = new LinearInterpolator();
            UnivariateFunction function = interpolator.interpolate(t0, P0);
            double[] P = new double[m];
            for (int i = 0; i < m; i++) {
                P[i] = function.value(t_seq[i]);
            }

            // initialize
            double[] X0 = pdf;

            // set up iteration
            double err = 1;
            int it = 0;
            double[] X = X0;

            // perform FFT
            double[] ft = padReal(pdf);
            FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
            Complex [] Ft = fft.transform(convertToComplex(ft), TransformType.FORWARD);

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
                    Xi[i] = X0[i] + 2 * I[i];
                }

                // compute error
                EuclideanDistance norm = new EuclideanDistance();
                err = norm.compute(X, Xi);

                // update
                X = Xi;
                it++;
            }

            if (it == maxit) {
                System.err.printf("calcB Warning: max iterations reached with error: %.2f%n", err);
            }

            // take last element
            B[x] = X[m - 1];
        }

        return B;
    }


    // Function for partial convolution using FFT
    protected static double[] convolveFFT(Complex[] fx, double[] y, int n, double eps) {

        // pad y
        double[] y_ext = padReal(y);

        // perform FFT on y_ext
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] fy = fft.transform(y_ext, TransformType.FORWARD);

        // element-wise multiplication of fx and fy (convolution in Fourier space)
        Complex[] fz = new Complex[fx.length];
        for (int i = 0; i < fx.length; i++) {
            fz[i] = fx[i].multiply(fy[i]);
        }

        // pad fz
        Complex [] fz_ext = padComplex(fz);

        // perform inverse FFT to get the result back in time domain
        Complex[] z = fft.transform(fz_ext, TransformType.INVERSE);

        // extract the real part and scale it by eps
        double[] z_real = new double[n];
        for (int i = 0; i < n; i++) {
            z_real[i] = z[i].getReal() * eps;
        }

        return z_real;
    }


    // Function for padding a real vector with 0 such that the length is a power of 2
    private static double[] padReal(double[] x) {
        int n = x.length;

        // find the next power of 2 greater than or equal to n
        int np = 1;
        while (np < n) {
            np *= 2;
        }

        // extend x
        double[] xp = new double[np];
        System.arraycopy(x, 0, xp, 0, n);

        return xp;
    }


    // Function for padding a complex vector with 0 such that the length is a power of 2
    private static Complex[] padComplex(Complex[] x) {
        int n = x.length;

        // find the next power of 2 greater than or equal to n
        int np = 1;
        while (np < n) {
            np *= 2;
        }

        // pad with 0
        Complex[] xp = new Complex[np];
        System.arraycopy(x, 0, xp, 0, n);
        for (int i = n; i < np; i++) {
            xp[i] = Complex.ZERO;
        }

        return xp;
    }
}
