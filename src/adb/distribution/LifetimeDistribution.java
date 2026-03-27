package adb.distribution;

import adb.util.Utils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;

import static adb.util.Utils.TRANSFORM_FORWARD;

public class LifetimeDistribution {
    // TODO: or implement as CalculationNode? (only needs to be re-calculated if lifetime or shape for a type change)

    double lifetime;
    double shape;
    protected Complex[] pdfFFT;
    protected double[] cdf;
    protected double[] seq;


    public LifetimeDistribution(double lifetime, double shape, double[] seq) {
        this.lifetime = lifetime;
        this.shape = shape;
        this.seq = seq;

        // initialize distribution
        double scale = lifetime / shape;
        GammaDistribution gammaDist = new GammaDistribution(shape, scale);

        // calculate PDF and CDF on given array
        int n = seq.length;
        double[] pdf = new double[n];
        cdf = new double[n];
        for (int i = 0; i < n; i++) {
            pdf[i] = Math.exp(gammaDist.logDensity(seq[i])); // use log to prevent underflow
            cdf[i] = gammaDist.cumulativeProbability(seq[i]);
        }

        // transform PDF
        pdfFFT = Utils.fft.transform(Utils.padZeros(pdf), TRANSFORM_FORWARD);
    }
}
