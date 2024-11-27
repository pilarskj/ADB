package test.adbp;

import adbp.GammaLogLikelihood;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

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
        double [] int_s = { 5, 8, 11 };
        double [] int_e = { 11, 11, 12 };
        double [] ext_e = { 5, 5, 8, 8 };

        // computational options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);

        double logL = GammaLogLikelihood.calcLogLikelihood(a, b, d, rho, origin, int_s, int_e, ext_e,
                maxIt, tolP, tolB, mP, mB, false);
        System.out.println(logL);

        assertEquals(logL, -26.91385,  0.05);
    }
}




