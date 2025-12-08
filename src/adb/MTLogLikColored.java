package adb;

import bdmmprime.distribution.SmallNumber;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import java.util.HashMap;
import java.util.stream.IntStream;

import static adb.GammaLogLikelihood.padZeros;
import static adb.MTLogLikelihood.calcMTP0;
import static org.apache.commons.math3.special.Gamma.logGamma;

public class MTLogLikColored {

    public static FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
    public static LinearInterpolator interpolator = new LinearInterpolator();

    public static double calcMTLogLikColored(double[] a, double[] b, double[] d, double rho,
                                             double[][] Xsi_s, double[][] Xsi_as,
                                             double t_or, int type_or, BranchList branches,
                                             int maxIt, double tolP, int mP) {

        // get number of types
        int ntypes = a.length;

        // generate linearly spaced values between 0 and origin
        int m = mP; // the number of time steps must be a power of 2 (required by FFT!)
        double[] seq = new double[m];
        double dx = t_or / m;
        for (int w = 0; w < m; w++) {
            seq[w] = dx * (w + 1);
        }
        assert seq[m - 1] == t_or;

        // calculate CDF and FFT from PDF from 0 to origin
        double[][] cdf = new double[m][ntypes];
        Complex[][] pdfFFT = new Complex[m*2][ntypes];
        IntStream.range(0, ntypes)
                .parallel()
                .forEach(i -> {
                    // get density and cumulative probability
                    double[] pdf = new double[m];
                    GammaDistribution gammaDist = new GammaDistribution(b[i], a[i]);
                    for (int w = 0; w < m; w++) {
                        pdf[w] = Math.exp(gammaDist.logDensity(seq[w]));
                        cdf[w][i] = gammaDist.cumulativeProbability(seq[w]);
                    }

                    // perform FFT
                    Complex[] Ft = fft.transform(padZeros(pdf), TransformType.FORWARD);
                    for (int w = 0; w < m*2; w++) {
                        pdfFFT[w][i] = Ft[w];
                    }
                });

        // calculate extinction probability over time
        double[][] P0 = calcMTP0(pdfFFT, cdf, seq, dx, d, rho, Xsi_s, Xsi_as, maxIt, tolP);

        // stop if extinction is certain
        if (P0[P0.length - 1][type_or] == 1.0) {
            return Double.NEGATIVE_INFINITY;
        }

        // extend tSeq to calculate probabilities of tiny branches
        double[] extSeq = new double[seq.length + 1];
        extSeq[0] = 0;
        System.arraycopy(seq, 0, extSeq, 1, seq.length);

        // extend and interpolate P0 (for calculating branch probabilities)
        HashMap<Integer, UnivariateFunction> P0Map = new HashMap<>();
        for (int i = 0; i < ntypes; i++) {
            double[] extP0 = new double[m + 1];
            extP0[0] = 1 - rho;
            for (int w = 0; w < m ; w++) {
                extP0[w + 1] = P0[w][i];
            }
            UnivariateFunction function = interpolator.interpolate(extSeq, extP0);
            P0Map.put(i, function);
        }

        // calculate probabilities of branch segments
        int nbranches = branches.listBranches().size();
        double[] branchProbs = new double[nbranches];
        for (int x = 0; x < nbranches; x++) {
            Branch branch = branches.getBranchByIndex(x);
            int type = branch.branchType;
            GammaDistribution gammaDist = new GammaDistribution(b[type], a[type]);
            if (branch.branchMode.equals("internal")) {
                branchProbs[x] = (1 - d[type]) * Math.exp(gammaDist.logDensity(branch.endTime - branch.startTime));
            } else {
                branchProbs[x] = rho * (1 - gammaDist.cumulativeProbability(branch.endTime));
            }
        }

        // start recursion
        SmallNumber treeL = calcLikelihoodBranchTyped(0, branches, Xsi_s, Xsi_as, ntypes, P0Map, branchProbs);

        // add tree factor
        int ntips = branches.countExternalBranches();
        double treeFactor = Math.log(2) * (ntips - 1) - logGamma(ntips + 1); // 2^(n-1)/n!

        return treeFactor - Math.log(1 - P0[P0.length - 1][type_or]) + treeL.log();
    }


    // Recursive function for calculating the tree likelihood (for branch-typed tree)
    public static SmallNumber calcLikelihoodBranchTyped(int branchIndex, BranchList branches, double[][] Xsi_s, double[][] Xsi_as, int ntypes,
                                                        HashMap<Integer, UnivariateFunction> P0Map, double[] branchProbs) {

        Branch branch = branches.getBranchByIndex(branchIndex);
        SmallNumber lik = new SmallNumber(branchProbs[branchIndex]);

        if (branch.branchMode.equals("external")) { // stop recursion at tips
            return lik;

        } else {
            int type = branch.branchType;

            if (branch.rightIndex > -1) { // two descendants
                int leftType = branches.getBranchByIndex(branch.leftIndex).branchType;
                int rightType = branches.getBranchByIndex(branch.rightIndex).branchType;
                SmallNumber leftLik = calcLikelihoodBranchTyped(branch.leftIndex, branches, Xsi_s, Xsi_as, ntypes, P0Map, branchProbs);
                SmallNumber rightLik = calcLikelihoodBranchTyped(branch.rightIndex, branches, Xsi_s, Xsi_as, ntypes, P0Map, branchProbs);
                double scalar;
                if (leftType == rightType) { // symmetric division
                    scalar = Xsi_s[type][leftType];
                } else { // asymmetric division
                    if (leftType == type) {
                        scalar = 0.5 * Xsi_as[type][rightType];
                    } else {
                        scalar = 0.5 * Xsi_as[type][leftType];
                    }
                }
                lik = lik.scalarMultiplyBy(scalar).multiplyBy(leftLik).multiplyBy(rightLik);

            } else { // one descendant
                double branchingTime = branch.startTime;
                int childType = branches.getBranchByIndex(branch.leftIndex).branchType;
                SmallNumber childLik = calcLikelihoodBranchTyped(branch.leftIndex, branches, Xsi_s, Xsi_as, ntypes, P0Map, branchProbs);
                double scalar; // TODO: factors 2 or 0.5?
                if (childType == type) { // i -> i: i -> ii or i -> ik
                    scalar = Xsi_s[type][type] * P0Map.get(type).value(branchingTime);
                    for (int i = 0; i < ntypes; i++) {
                        scalar += Xsi_as[type][i] * P0Map.get(i).value(branchingTime);
                    }
                } else { // i -> j: i -> jj or i -> ij
                    scalar = Xsi_s[type][childType] * P0Map.get(childType).value(branchingTime) +
                            Xsi_as[type][childType] * P0Map.get(type).value(branchingTime);
                }
                lik = lik.scalarMultiplyBy(scalar).multiplyBy(childLik);
            }
            return lik;
        }
    }

}
