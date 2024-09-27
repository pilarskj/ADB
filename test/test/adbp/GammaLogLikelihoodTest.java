package test.adbp;

import adbp.GammaLogLikelihood;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class GammaLogLikelihoodTest {

    @Test
    public void testFFT() {

        // Create a FastFourierTransformer instance
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);

        // Define an easy input: a simple array of real numbers
        double[] input = {1.0, 2.0, 3.0, 4.0};

        // Perform the forward FFT transform
        Complex[] fftResult = fft.transform(input, TransformType.FORWARD);

        // Print the FFT result
        System.out.println("FFT result:");
        for (Complex c : fftResult) {
            System.out.println(c);
        }

        // Perform the inverse FFT transform to get back the original data
        Complex[] inverseFftResult = fft.transform(fftResult, TransformType.INVERSE);

        // Print the inverse FFT result (it should closely match the original input)
        System.out.println("\nInverse FFT result:");
        for (Complex c : inverseFftResult) {
            System.out.println(c);
        }
    }

    @Test
    public void testConvolveFFT() {

        // Example input for fx (complex vector)
        Complex[] fx = {
                new Complex(1, 2),
                new Complex(3, 4),
                new Complex(5, 6)
        };

        // Example input for y (real vector)
        double[] y = { 1.0, 2.0, 3.0, 4.0, 4.0 };

        // Parameters for the convolution
        int n = 4;
        double eps = 0.1;  // Scaling factor

        // Perform convolution using FFT
        double[] result = GammaLogLikelihood.convolveFFT(fx, y, n, eps);

        // Output the result
        for (double value : result) {
            System.out.println(value);
        }
    }

    @Test
    public void testP0() {

        // Example input
        double rho = 0.6;
        double a = 0.03;
        double b = 570;
        double t_or = 100;

        int m = 10;
        double[] t = new double[m];
        double dx = t_or / (m - 1);
        for (int i = 0; i < m; i++) {
            t[i] = dx * i;
        }
        assert t[m - 1] == t_or;

        int maxit = 100;

        double[] result = GammaLogLikelihood.calcP0(rho, a, b, t, dx, maxit);

        for (double value : result) {
            System.out.println(value);
        }
    }

    @Test
    public void testLogLikelihood() {

        // cf. R:
        // tree = read.tree(text = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;")
        // times = get_tree_times(tree = tree, t_or = 11) # does not include the stem branch!
        // (extract branches)

        // inputs
        // phylodynamic parameters
        double rho = 0.1;
        double a = 1;
        double b = 1;
        double t_or = 12;
        // branch lengths
        double [] int_s = { 5, 8, 11 };
        double [] int_e = { 11, 11, 12 };
        double [] ext_e = { 5, 5, 8, 8 };
        // computational options
        int m = (int)Math.pow(2, 14);
        int maxit = 100;

        double logL = GammaLogLikelihood.calcLogLikelihood(rho, a, b, t_or, int_s, int_e, ext_e, m, maxit);
        System.out.println(logL);

        // assertEquals(logL, -25.8506,  0.05);
    }
}




