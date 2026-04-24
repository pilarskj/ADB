package adb.distribution;

import adb.util.Utils;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;

import java.util.BitSet;

import static adb.util.Utils.TRANSFORM_FORWARD;


public class LifetimeDistributions extends CalculationNode {
    // TODO: test storing/ restoring in MCMC

    public Input<RealParameter> lifetimeParameterInput =
            new Input<>("lifetime", "", Input.Validate.REQUIRED);
    public Input<Parameter> shapeParameterInput =
            new Input<>("shape", "", Input.Validate.REQUIRED);

    public Input<double[]> timeArrayInput =
             new Input<>("timeArray", "", Input.Validate.REQUIRED);

    RealParameter lifetimeParameter;
    Parameter shapeParameter;
    int nTypes;
    double[] timeArray;

    // for identifying dirty indices
    boolean dirty;
    BitSet dirtyIndices = new BitSet();

    // output
    LifetimeDistribution[] lifetimeDistributions;
    LifetimeDistribution[] storedLifetimeDistributions;


    // Internal class
    public static class LifetimeDistribution {
        double lifetime;
        double shape;
        Complex[] pdfFFT;
        double[] cdf;

        public LifetimeDistribution(double lifetime, double shape, Complex[] pdfFFT, double[] cdf) {
            this.lifetime = lifetime;
            this.shape = shape;
            this.pdfFFT = pdfFFT;
            this.cdf = cdf;
        }

        public void copyFrom(LifetimeDistribution other) {
            // deep copy of all parameters
            this.lifetime = other.lifetime;
            this.shape = other.shape;
            for (int i = 0; i < other.pdfFFT.length; i++) {
                this.pdfFFT[i] = other.pdfFFT[i];
            }
            System.arraycopy(other.cdf, 0, this.cdf, 0, other.cdf.length);
        }

        public double[] getCDF() {
            return cdf;
        }

        public Complex[] getTransformedPDF() {
            return pdfFFT;
        }
    }


    @Override
    public void initAndValidate() {
        lifetimeParameter = lifetimeParameterInput.get();
        shapeParameter = shapeParameterInput.get();
        nTypes = lifetimeParameter.getDimension();
        timeArray = timeArrayInput.get();

        lifetimeDistributions = new LifetimeDistribution[nTypes];
        storedLifetimeDistributions = new LifetimeDistribution[nTypes];

        // get all initial distributions
        dirty = true;
        dirtyIndices.set(0, nTypes);
        update();
        dirtyIndices.clear();

        // copy to stored (avoid null)
        for (int i = 0; i < nTypes; i++) {
            storedLifetimeDistributions[i] = new LifetimeDistribution(
                    lifetimeDistributions[i].lifetime,
                    lifetimeDistributions[i].shape,
                    lifetimeDistributions[i].pdfFFT,
                    lifetimeDistributions[i].cdf);
        }
    }


    private void update() {
        if (!dirty) return;

        // only update what actually changed!
        for (int i = dirtyIndices.nextSetBit(0); i >= 0; i = dirtyIndices.nextSetBit(i+1)) {
            updateLifetimeDistribution(i);
        }

        dirty = false;
    }


    private void updateLifetimeDistribution(int x) {
        double lifetime = lifetimeParameter.getArrayValue(x);
        double shape = shapeParameter.getArrayValue(x);

        // initialize distribution
        double scale = lifetime / shape;
        GammaDistribution gammaDist = new GammaDistribution(shape, scale);

        // calculate PDF and CDF on given array
        int n = timeArray.length;
        double[] pdf = new double[n];
        double[] cdf = new double[n];
        for (int i = 0; i < n; i++) {
            pdf[i] = Math.exp(gammaDist.logDensity(timeArray[i])); // use log to prevent underflow
            cdf[i] = gammaDist.cumulativeProbability(timeArray[i]);
        }

        // transform PDF
        Complex[] pdfFFT = Utils.fft.transform(Utils.padZeros(pdf), TRANSFORM_FORWARD);

        // store in list
        lifetimeDistributions[x] = new LifetimeDistribution(lifetime, shape, pdfFFT, cdf);
    }


    public LifetimeDistribution[] getLifetimeDistributions() {
        update();
        return lifetimeDistributions;
    }


    @Override
    protected boolean requiresRecalculation() {
        for (int i = 0; i < nTypes; i++) {
            if (lifetimeParameter.isDirty(i)) { // Note: for some reason shapeParameter.isDirty(i) does not work!
                dirtyIndices.set(i); // mark this index as needing a store/restore
            }
            if (shapeParameter instanceof RealParameter) {
                if (((RealParameter)shapeParameter).isDirty(i)) {
                    dirtyIndices.set(i);
                }
            } else if (shapeParameter instanceof IntegerParameter) {
                if (((IntegerParameter)shapeParameter).isDirty(i)) {
                    dirtyIndices.set(i);
                }
            }
        }
        dirty = true;
        return true;
    }


    @Override
    protected void store() {
        // only iterate over bits that are set to true
        for (int i = dirtyIndices.nextSetBit(0); i >= 0; i = dirtyIndices.nextSetBit(i+1)) {
            storedLifetimeDistributions[i].copyFrom(lifetimeDistributions[i]);
        }
        super.store();
    }


    @Override
    protected void restore() {
        for (int i = dirtyIndices.nextSetBit(0); i >= 0; i = dirtyIndices.nextSetBit(i+1)) {
            lifetimeDistributions[i].copyFrom(storedLifetimeDistributions[i]);
        }
        dirtyIndices.clear();
        super.restore();
    }


    // TODO: necessary? -- check clearing of dirtyIndices
    @Override
    protected void accept() {
        /* for (int i = dirtyIndices.nextSetBit(0); i >= 0; i = dirtyIndices.nextSetBit(i+1)) {
            storedLifetimes[i] = lifetimeParameter.getArrayValue(i);
            storedShapes[i] = shapeParameter.getArrayValue(i);
            storedLifetimeDistributions[i].copyFrom(lifetimeDistributions[i]);
        } */
        dirtyIndices.clear();
        super.accept();
    }
}
