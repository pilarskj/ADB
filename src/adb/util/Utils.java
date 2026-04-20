package adb.util;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class Utils {

    // Global tools for calculations
    public static FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
    public static LinearInterpolator interpolator = new LinearInterpolator();

    public static TransformType TRANSFORM_FORWARD = TransformType.FORWARD;
    public static TransformType TRANSFORM_INVERSE = TransformType.INVERSE;


    // Check if number is a power of 2
    // https://stackoverflow.com/questions/600293/how-to-check-if-a-number-is-a-power-of-2
    public static boolean isPowerOfTwo(int x) {
        return (x > 0) && ((x & (x - 1)) == 0);
    }


    // Pad an array with 0 to its double length
    public static double[] padZeros(double[] x) {
        int n = x.length;
        double[] xp = new double[n * 2];
        System.arraycopy(x, 0, xp, 0, n);
        return xp;
    }


    // Generate linearly spaced array (excluding start, including end)
    public static double[] linSpace(double start, double end, int n) {
        double dx = (end - start) / n;
        double[] seq = new double[n];
        for (int i = 0; i < n; i++) {
            seq[i] = start + dx * (i + 1);
        }
        return seq;
    }


    // Perform partial convolution using FFT
    public static double[] convolveFFT(Complex[] fx, double[] y, int n, double eps) {

        // perform FFT on padded y
        Complex[] fy = fft.transform(padZeros(y), TRANSFORM_FORWARD);

        // element-wise multiplication of fx and fy (convolution in Fourier space)
        Complex[] fz = new Complex[fx.length];
        for (int i = 0; i < fx.length; i++) {
            fz[i] = fx[i].multiply(fy[i]);
        }

        // perform inverse FFT to get the result back in time domain
        Complex[] z = fft.transform(fz, TRANSFORM_INVERSE);

        // extract the real part and scale it by eps
        double[] z_real = new double[n];
        for (int i = 0; i < n; i++) {
            z_real[i] = z[i].getReal() * eps;
        }

        return z_real;
    }


    // Calculate the mean of an array
    public static double getMean(double[] array) {
        double sum = 0;
        for (double num : array) {
            sum += num;
        }
        return sum / array.length;
    }


    // Enforce non-increasing or non-decreasing array -- change String order?
    public static void forceOrder(double[] array, String order) {
        if (order.equals("non-increasing")) {
            for (int i = 1; i < array.length; i++) {
                if (array[i] > array[i - 1]) {
                    array[i] = array[i - 1];  // adjust to maintain non-increasing order
                }
            }
        } else if (order.equals("non-decreasing")) {
            for (int i = 1; i < array.length; i++) {
                if (array[i] < array[i - 1]) {
                    array[i] = array[i - 1];  // adjust to maintain non-decreasing order
                }
            }
        } else {
            throw new IllegalArgumentException("Invalid order: use non-increasing or non-decreasing");
        }
    }


    // Find the index of the closest value in a sorted array using binary search // TODO: overcomplicated?
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


    // Calculate L2 distance element-wise
    public static double getError(double[] x, double[] y) {
        int n = x.length;
        double sum = 0;
        for (int i = 0; i < n; i++) {
            double diff = x[i] - y[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    public static double getError(double[][] X, double[][] Y) {
        int n = X.length; // number of rows
        int m = X[0].length;  // number of columns
        double sum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double diff = X[i][j] - Y[i][j];
                sum += diff * diff;
            }
        }
        return Math.sqrt(sum);
    }

    public static double getError(double[][][] X, double[][][] Y){
        int n = X.length; // number of arrays
        int m = X[0].length;  // number of rows
        int o = X[0][0].length;  // number of columns
        double sum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                for (int k = 0; k < o; k++) {
                    double diff = X[i][j][k] - Y[i][j][k];
                    sum += diff * diff;
                }
            }
        }
        return Math.sqrt(sum);
    }


    public static void saveArrays(double[] x, double[] y, int thin, String header, String fileName) {
        try (PrintWriter writer = new PrintWriter(new FileWriter(fileName))) {
            writer.println(header);

            for (int i = 0; i < x.length; i += thin) {
                writer.printf("%.6f,%.6f%n", x[i], y[i]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
