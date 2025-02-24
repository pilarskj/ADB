package test.adbp;

import adbp.GammaLogLikelihood;
import beast.base.inference.parameter.IntegerParameter;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Arrays;

import static adbp.GammaLogLikelihood.findClosestIndex;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class GammaLogLikelihoodTest {

    @Test
    public void testGammaLogLikelihood() {

        // phylodynamic parameters
        double a = 1;
        int b = 1;
        double d = 0;
        double rho = 0.1;
        double origin = 12;

        // branch lengths for tree = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;"
        double [] intS = { 5, 8, 11 };
        double [] intE = { 11, 11, 12 };
        double [] extE = { 5, 5, 8, 8 };

        // computational options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);

        double logL = GammaLogLikelihood.calcLogLikelihood(a, b, d, rho, origin, intS, intE, extE,
                maxIt, tolP, tolB, mP, mB, false);
        System.out.println(logL);

        assertEquals(logL, -26.91385,  0.05);
    }

    @Test
    public void testFindClosestIndex() {

        double[] array = {1, 1.5, 3, 4.5, 6, 7.5, 9};
        double x = 0.5;
        double y = 10;

        int ix = findClosestIndex(array, x);
        int iy = findClosestIndex(array, y);
        System.out.println(ix + " " + iy);
    }

    @Test
    public void testP0() throws Exception {

        // phylodynamic parameters
        double a = 0.05;
        Number b = 100;
        double d = 0.1;
        double rho = 0.0001;
        double origin = 100;

        // branch lengths for tree = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;"
        double [] intS = { 5, 8, 11 };
        double [] intE = { 11, 11, 12 };
        double [] extE = { 5, 5, 8, 8 };

        // computational options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);

        double logL = GammaLogLikelihood.calcLogLikelihood(a, b, d, rho, origin, intS, intE, extE,
                maxIt, tolP, tolB, mP, mB, false);
    }
}




