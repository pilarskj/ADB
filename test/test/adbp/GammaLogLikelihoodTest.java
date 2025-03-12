package test.adbp;

import adbp.GammaLogLikelihood;
import org.junit.jupiter.api.Test;

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

        // analytical BD estimate
        double bd = GammaLogLikelihood.calcLogLikelihood(a, b, d, rho, origin, intS, intE, extE,
                maxIt, tolP, tolB, mP, mB, false, true);
        System.out.println("BD logL = " + bd);

        // ADB exact FFT solution
        double adb = GammaLogLikelihood.calcLogLikelihood(a, b, d, rho, origin, intS, intE, extE,
                maxIt, tolP, tolB, mP, mB, false, false);
        System.out.println("ADB logL = " + adb);

        assertEquals(adb, bd,  0.05);
    }
}




