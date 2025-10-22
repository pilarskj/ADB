package adb;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/* Adapted from beast.base.inference.operator.DeltaExchangeOperator */
@Description("Operator for row-wise sum-constrained parameters.")
public class MatrixDeltaExchangeOperator extends Operator {

    public Input<List<RealParameter>> parameterInput =
            new Input<>("parameter", "One or more parameters to operate on", new ArrayList<>());

    public Input<Double> deltaInput =
            new Input<>("delta", "Magnitude of change for two randomly picked values.", 0.1);

    public Input<Boolean> excludeZerosInput =
            new Input<>("excludeZeros", "Exclude zeros from operation (default true)", true);

    public final Input<Boolean> autoOptimizeInput =
            new Input<>("autoOptimize", "Adjust delta during the MCMC run to improve mixing (default true).", true);


    private List<RealParameter> parameters;
    private double delta;
    private boolean excludeZeros;
    private boolean autoOptimize;

    private int nparams;
    private int ndims;
    private List<double[][]> matrices;
    private List<Map<Integer, int[]>> indices;


    public void initAndValidate() {
        parameters = parameterInput.get();
        delta = deltaInput.get();
        excludeZeros = excludeZerosInput.get();
        autoOptimize = autoOptimizeInput.get();

        nparams = parameters.size();
        int dim = parameters.get(0).getDimension();
        ndims = (int)(Math.sqrt(dim));

        // if more parameters are provided, assert that they have the same dimension
        if (nparams > 1) {
            for (RealParameter param : parameters) {
                if (param.getDimension() != dim) {
                    throw new IllegalArgumentException("Parameter dimensions do not match.");
                }
            }
        }

        // collect non-zero elements
        if (excludeZeros) { // TODO: Am I overcomplicating?
            matrices = new ArrayList<>(nparams);
            indices = new ArrayList<>(nparams);
            List<int[]> nonZeroCounts = new ArrayList<>(nparams);

            for (RealParameter param : parameters) {
                // populate matrices
                double[][] matrix = new double[ndims][ndims];
                for (int i = 0; i < ndims; i++) {
                    for (int j = 0; j < ndims; j++) {
                        matrix[i][j] = param.getDoubleValues()[i * ndims + j];
                    }
                }
                matrices.add(matrix);

                // collect 0 indices row-wise
                int[] counts = new int[ndims];
                Map<Integer, int[]> idx = new HashMap<>();
                for (int i = 0; i < ndims; i++) {

                    // first, count non-zero elements
                    int count = 0;
                    for (double v : matrix[i]) if (v != 0.0) count++;
                    counts[i] = count;

                    // store indices
                    if (count > 0) {
                        int[] ix = new int[count];
                        int k = 0;
                        for (int j = 0; j < ndims; j++) {
                            if (matrix[i][j] != 0.0) ix[k++] = j;
                        }
                        idx.put(i, ix);
                    }
                }
                indices.add(idx);
                nonZeroCounts.add(counts);
            }

            // check valid rows
            boolean valid = false;
            for (int i = 0; i < ndims; i++) {
                int rowCount = 0;
                for (int p = 0; p < nparams; p++) {
                    rowCount += nonZeroCounts.get(p)[i];
                }
                if (rowCount > 1) { valid = true; break; } // at least one row can be modified with this operator
            }
            if (!valid) {
                throw new IllegalArgumentException("Provided matrices are too sparse for this operator.");
            }
        }
    }


    @Override
    public final double proposal() {

        // select matrices and entries to modify
        int m1 = Randomizer.nextInt(nparams);
        int m2 = Randomizer.nextInt(nparams);
        int i = Randomizer.nextInt(ndims); // row

        // columns
        int j1;
        int j2;
        if (!excludeZeros) {
            j1 = Randomizer.nextInt(ndims);
            j2 = Randomizer.nextInt(ndims);
        } else {
            int k = i;
            while (indices.get(m1).get(i) == null || indices.get(m2).get(i) == null) { //
                if (i < ndims-1) {i++;} else {i=0;}
                if (i == k) { return Double.NEGATIVE_INFINITY; } // no move possible; or 0.0? // TODO: assure correctness for all cases!
            }
            j1 = randomSample(indices.get(m1).get(i));
            j2 = randomSample(indices.get(m2).get(i));
        }

        int dim1 = i * ndims + j1; // entry 1
        int dim2 = i * ndims + j2; // entry 2
        if (m1 == m2 && dim1 == dim2) { // no move
            return 0.0;
        }

        // extract values
        double value1 = parameters.get(m1).getValue(dim1);
        double value2 = parameters.get(m2).getValue(dim2);

        // update values
        double d = Randomizer.nextDouble() * delta; // d ~ Uniform(0, delta)
        value1 -= d / 2;
        value2 += d / 2;

        if (value1 < parameters.get(m1).getLower() || value1 > parameters.get(m1).getUpper() ||
                value2 < parameters.get(m2).getLower() || value2 > parameters.get(m2).getUpper()) {
            return Double.NEGATIVE_INFINITY;
        } else {
            parameters.get(m1).setValue(dim1, value1);
            parameters.get(m2).setValue(dim2, value2);
        }

        return 0.0; // symmetric move
    }


    public static int randomSample(int[] values) {
        int idx = Randomizer.nextInt(values.length);
        return values[idx];
    }


    @Override
    public double getCoercableParameterValue() {
        return delta;
    }


    @Override
    public void setCoercableParameterValue(double value) {
        delta = value;
    }


    @Override
    public void optimize(double logAlpha) {
        if (autoOptimize) {
            double _delta = calcDelta(logAlpha);
            _delta += Math.log(delta);
            delta = Math.exp(_delta);
        }
    }


    @Override
    public final String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new magnitude of change
        double newDelta = delta * ratio;
        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else return "";
    }


    @Override
    public List<StateNode> listStateNodes() {
        return new ArrayList<>(parameterInput.get());
    }
}
