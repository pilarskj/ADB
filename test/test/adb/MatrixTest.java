package test.adb;

import beast.base.util.Randomizer;
import org.junit.jupiter.api.Test;

import java.util.*;

import static adb.MatrixDeltaExchangeOperator.randomSample;

public class MatrixTest {

    @Test
    public void testMatrix() {

        int nparams = 2;
        int ndims = 3;

        List<double[][]> matrices = new ArrayList<>(nparams);
        double[][] matrix1 = {{1,1,0}, {0,0,0},{1,0,1}};
        //double[][] matrix1 = {{0.5,0},{0,0.5}};
        matrices.add(matrix1);
        double[][] matrix2 = {{0,0,0}, {1,0,0}, {0,1,0}};
        //double[][] matrix2 = {{0,0.5},{0.5,0}};
        matrices.add(matrix2);

        List<int[]> nonZeroCounts = new ArrayList<>(nparams);
        List<Map<Integer, int[]>> indices = new ArrayList<>(nparams);
        for (int m = 0; m < nparams; m++) {
            double[][] matrix = matrices.get(m);

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
                // System.out.println("Row " + i + ": " + Arrays.toString(ix));
            }
            indices.add(idx);
            nonZeroCounts.add(counts);
        }

        // collect valid rows
        boolean valid = false;
        boolean[] validRows = new boolean[ndims];
        for (int i = 0; i < ndims; i++) {
            int rowCount = 0;
            for (int p = 0; p < nparams; p++) {
                rowCount += nonZeroCounts.get(p)[i];
            }
            if (rowCount > 1) {
                valid = true;
                validRows[i] = true;
            } else { validRows[i] = false; }
        }


        // select matrices and entries to modify
        int m1 = Randomizer.nextInt(nparams);
        int m2 = Randomizer.nextInt(nparams);
        int i = Randomizer.nextInt(ndims); // row
        int k = i;
        while (indices.get(m1).get(i) == null || indices.get(m2).get(i) == null) {
            if (i < ndims-1) {i++;} else {i=0;}
            if (i == k) { System.out.println("No move possible!"); } // no move possible; or Double.NEGATIVE_INFINITY?
        }
        int j1 = randomSample(indices.get(m1).get(i));
        int j2 = randomSample(indices.get(m2).get(i));

        int dim1 = i * ndims + j1; // entry 1
        int dim2 = i * ndims + j2; // entry 2

        System.out.println("Matrix " + m1 + ", entry " + dim1);
        System.out.println("Matrix " + m2 + ", entry " + dim2);
    }
}
