import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.junit.jupiter.api.Test;

public class FFTTest {

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
        double[] result = GammaBranchingModel.convolveFFT(fx, y, n, eps);

        // Output the result
        for (double value : result) {
            System.out.println(value);
        }
    }

    @Test
    public void testCalcP0() {

        // Example input
        double rho = 0.5;
        double a = 2;
        double b = 1;
        double[] t = { 6, 4, 5, 8.2, 4.1, 9.8 };
        double dx = 1;
        int maxit = 5;

        double[] result = GammaBranchingModel.calcP0(rho, a, b, t, dx, maxit);

        for (double value : result) {
            System.out.println(value);
        }
    }
}




