import beast.base.evolution.tree.TreeDistribution;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.UnivariateFunction;

import static org.apache.commons.math3.complex.ComplexUtils.convertToComplex;

public class GammaBranchingModel {// extends TreeDistribution {

    // Function for calculating the probability of a single descendant
    protected static double[] calcP1(double rho, double a, double b, double[] t, double[] t0, double[] P0, double dx, int maxit) {

        // get length
        int n = t0.length;

        // calculate the PDF and CDF of the gamma distribution
        double[] pdf = new double[n];
        double[] cdf = new double[n];
        GammaDistribution gammaDist = new GammaDistribution(b, a);
        for (int i = 0; i < n; i++) {
            pdf[i] = gammaDist.density(t[i]);
            cdf[i] = gammaDist.cumulativeProbability(t[i]);
        }

        // initialize L
        double[] L0 = new double[n];
        for (int i = 0; i < n; i++) {
            L0[i] = rho * (1.0 - cdf[i]);
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[] L = L0;

        // perform FFT
        double[] ft = padReal(pdf);
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex [] Ft = fft.transform(convertToComplex(ft), TransformType.FORWARD);

        // iterate
        while (err > 1e-12 && it < maxit) {

            // multiply
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = P0[i] * L[i];
            }

            // partially convolve
            double[] I = convolveFFT(Ft, y, n, dx);

            // sum
            double[] Li = new double[n];
            for (int i = 0; i < n; i++) {
                Li[i] = L0[i] + 2 * I[i];
            }

            // compute error
            EuclideanDistance norm = new EuclideanDistance();
            err = norm.compute(L, Li);

            // update
            L = Li;
            it++;
        }
        // System.err.printf("Warning: max iterations reached with error: %.2f%n", err);

        // interpolate
        LinearInterpolator interpolator = new LinearInterpolator();
        UnivariateFunction function = interpolator.interpolate(t0, L);
        double[] P1 = new double[t.length];
        for (int i = 0; i < t.length; i++) {
            P1[i] = function.value(t[i]);
        }

        return P1;
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

        // initialize P00
        double[] P00 = new double[n];
        for (int i = 0; i < n; i++) {
            P00[i] = (1.0 - rho) * (1.0 - cdf[i]);
        }

        // set up iteration
        double err = 1;
        int it = 0;
        double[] P0 = P00;

        // perform FFT
        double[] ft = padReal(pdf);
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] Ft = fft.transform(convertToComplex(ft), TransformType.FORWARD);

        // iterate
        while (err > 1e-12 && it < maxit) {

            // square
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                y[i] = P0[i] * P0[i];
            }

            // partially convolve
            double[] I = convolveFFT(Ft, y, n, dx);

            // sum
            double[] Pi = new double[n];
            for (int i = 0; i < n; i++) {
                Pi[i] = P00[i] + I[i];
            }

            // compute error
            EuclideanDistance norm = new EuclideanDistance();
            err = norm.compute(Pi, P0);

            // update
            P0 = Pi;
            it++;
        }
        // System.err.printf("Warning: max iterations reached with error: %.2f%n", err);
        return P0;
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
