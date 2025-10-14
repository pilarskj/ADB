package adb;

import bdmmprime.distribution.SmallNumber;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.Pair;

import java.util.HashMap;
import java.util.stream.IntStream;

import static adb.GammaLogLikelihood.convolveFFT;
import static adb.GammaLogLikelihood.padZeros;
import static org.apache.commons.math3.special.Gamma.logGamma;

/*
Class for solving the equations for P0, P1 and B in the multi-type case
and for calculating the likelihood of a tip-typed tree
(based on the parameters, branching times, and types at tips).
 */
public class MTLogLikelihood {

    public static FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
    public static LinearInterpolator interpolator = new LinearInterpolator();

    public static double calcMTLogLikelihood(double[] a, double[] b, double[] d, double rho,
                                             double[][] Xsi_s, double[][] Xsi_as,
                                             double t_or, int type_or, BranchList branches,
                                             int maxIt, double tolP, double tolB, int mP, int mB) {

        // get number of types
        int ntypes = a.length;
        assert b.length == ntypes; // TODO: put assertions into SpeciesTreeDistribution class
        assert d.length == ntypes;

        // assert type at origin
        assert branches.getBranchByIndex(0).endType == type_or;

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

        /* // print (subsequence)
        for (int w = 0; w < m; w += 128) {
            System.out.println(seq[w] + "," + P0[w][0] + "," + P0[w][1]);
        } */

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

        /* // evaluate P0Map at branching times
        for (Branch branch : branches.listBranches()) {
            System.out.println(branch.endTime + "," + P0Map.get(0).value(branch.endTime) + "," + P0Map.get(1).value(branch.endTime) + ",ADB");
        } */

        // calculate probability of single descendant over time
        double[][][] P1 = calcMTP1(pdfFFT, cdf, seq, dx, P0, d, rho, Xsi_s, Xsi_as, maxIt, tolP);

        // interpolate P1
        HashMap<Pair<Integer,Integer>, UnivariateFunction> P1Map = new HashMap<>();
        for (int i = 0; i < ntypes; i++) {
            for (int j = 0; j < ntypes; j++) {
                // extract array from P1 (and extend to 0)
                double[] extP1 = new double[m + 1];
                if (i == j) { extP1[0] = rho; } else { extP1[0] = 0; }
                for (int w = 0; w < m ; w++) {
                    extP1[w + 1] = P1[w][i][j];
                }
                // interpolate
                UnivariateFunction function = interpolator.interpolate(extSeq, extP1);
                P1Map.put(new Pair<>(i, j), function);
            }
        }

        /* // print functions over time array
        for (int w = 0; w < m; w += 128) {
            System.out.println(seq[w] + "," +
                    P1Map.get(new Pair<>(0,0)).value(seq[w]) + "," +
                    P1Map.get(new Pair<>(0,1)).value(seq[w]) + "," +
                    P1Map.get(new Pair<>(1,0)).value(seq[w]) + "," +
                    P1Map.get(new Pair<>(1,1)).value(seq[w]));
        } */

        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();

        // collect times of internal branches
        double[] intS = new double[n_int];
        double[] intE = new double[n_int];
        for (int x = 0; x < n_int; x++) {
            Branch branch = branches.getBranchByIndex(x);
            intS[x] = branch.startTime;
            intE[x] = branch.endTime;
        }

        // calculate probabilities of internal branches
        double[][][] B = calcMTB(a, b, d, Xsi_s, Xsi_as, intS, intE, P0Map, maxIt, tolB, mB);

        // store subtree-likelihoods (dynamic programming)
        SmallNumber[][] subtreeLik = new SmallNumber[n_int + n_ext][ntypes];

        // start recursion
        SmallNumber treeL = calcLikelihoodTipTyped(branches, n_ext, ntypes, 0, type_or,
                Xsi_s, Xsi_as, P1Map, B, subtreeLik);

        /* // print output (debugging): node, ge0, ge1
        for (int x = 0; x < subtreeLik.length; x++) {
            Branch branch = branches.getBranchByIndex(x);
            if (subtreeLik[x][1] == null) { subtreeLik[x][1] = new SmallNumber(0.0); }
            System.out.println(x + "," + branch.startNode + "," + subtreeLik[x][0].toString() + "," + subtreeLik[x][1].toString()
                    + "," + branch.endTime + "," + branch.startTime + "," + branch.branchMode);
        } */

        // add tree factor
        double treeFactor = Math.log(2) * (n_ext - 1) - logGamma(n_ext + 1); // 2^(n-1)/n!

        return treeFactor - Math.log(1 - P0[P0.length - 1][type_or]) + treeL.log();
    }


    // Recursive function for calculating the tree likelihood (for tip-typed tree)
    public static SmallNumber calcLikelihoodTipTyped(BranchList branches, int ntips, int ntypes,
                                                     int stem_index, int stem_type, double[][] Xsi_s, double[][] Xsi_as,
                                                     HashMap<Pair<Integer,Integer>, UnivariateFunction> P1Map, double[][][] B,
                                                     SmallNumber[][] subtreeLik) {

        // check if calculation is stored already
        if (subtreeLik[stem_index][stem_type] != null) {
            return subtreeLik[stem_index][stem_type];
        }

        SmallNumber result;

        Branch branch = branches.getBranchByIndex(stem_index);
        // at tips
        if (stem_index >= ntips - 1) {
            result = new SmallNumber(P1Map.get(new Pair<>(stem_type, branch.startType)).value(branch.endTime));

        } else {
            // recursion
            SmallNumber subtreeL = new SmallNumber();

            for (int j = 0; j < ntypes; j++){
                SmallNumber branchL = new SmallNumber(B[stem_index][stem_type][j]);

                // get subtrees
                int right_stem = branch.rightIndex;
                int left_stem = branch.leftIndex;
                SmallNumber right_j = calcLikelihoodTipTyped(branches, ntips, ntypes, right_stem, j, Xsi_s, Xsi_as, P1Map, B, subtreeLik);
                SmallNumber left_j = calcLikelihoodTipTyped(branches, ntips, ntypes, left_stem, j, Xsi_s, Xsi_as, P1Map, B, subtreeLik);

                // loop over possible type transitions at internal nodes
                SmallNumber sumL = new SmallNumber();
                for (int k = 0; k < ntypes; k++){
                    SmallNumber right_k = calcLikelihoodTipTyped(branches, ntips, ntypes, right_stem, k, Xsi_s, Xsi_as, P1Map, B, subtreeLik);
                    SmallNumber left_k = calcLikelihoodTipTyped(branches, ntips, ntypes, left_stem, k, Xsi_s, Xsi_as, P1Map, B, subtreeLik);

                    sumL = sumL.addTo((left_k.multiplyBy(right_k)).scalarMultiplyBy(Xsi_s[j][k]));
                    SmallNumber optionA = left_k.multiplyBy(right_j);
                    SmallNumber optionB = left_j.multiplyBy(right_k);
                    sumL = sumL.addTo((optionA.addTo(optionB)).scalarMultiplyBy(0.5 * Xsi_as[j][k]));
                }

                subtreeL = subtreeL.addTo(branchL.multiplyBy(sumL));
            }
            result = subtreeL;
        }
        subtreeLik[stem_index][stem_type] = result;
        return result;
    }


    // Function for calculating the extinction probability
    public static double[][] calcMTP0(Complex[][] pdfFFT, double[][] cdf, double[] seq, double dx,
                                      double[] d, double rho, double[][] Xsi_s, double[][] Xsi_as,
                                      int maxIt, double tol) {

        // notation: it = iteration, w = integration variable (time), i,j,k = types
        // get number of types and time steps
        int ntypes = d.length;
        int m = seq.length;

        // initialize matrix
        double[][] X0 = new double[m][ntypes];
        IntStream.range(0, ntypes)
                .parallel()
                .forEach(i -> {
                    for (int w = 0; w < m; w++) {
                        X0[w][i] = (1 - rho) * (1 - cdf[w][i]) + d[i] * cdf[w][i];
                    }
                });

        // set up iteration
        double err = 1;
        int it = 0;
        double[][] X = X0;

        // iterate
        while (err > tol && it < maxIt) {
            double[][] Xi = new double[m][ntypes];

            for (int i = 0; i < ntypes; i++) {
                // get vectors for convolution
                double[] y = new double[m];
                for (int w = 0; w < m; w++) { // multiply elementwise on times
                    for (int j = 0; j < ntypes; j++) { // sum over all types k
                        y[w] += Xsi_s[i][j] * X[w][j] * X[w][j] + Xsi_as[i][j] * X[w][i] * X[w][j];
                    }
                }
                // extract column from the pdf matrix
                Complex[] Ft = new Complex[m*2];
                for (int w = 0; w < m*2; w++) {
                    Ft[w] = pdfFFT[w][i];
                }

                // partially convolve
                double[] I = convolveFFT(Ft, y, m, dx);

                // sum
                for (int w = 0; w < m; w++) {
                    Xi[w][i] = X0[w][i] + (1 - d[i]) * I[w];
                }
            }

            // compute error
            err = getMatrixError(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxIt) {
            System.err.printf("calcP0 Warning: max iterations reached with error: %.2f%n", err);
        }

        return X;
    }


    // Function for calculating the probability of a single descendant (considering start and end type)
    public static double[][][] calcMTP1(Complex[][] pdfFFT, double[][] cdf, double[] seq, double dx, double[][] P0,
                                        double[] d, double rho, double[][] Xsi_s, double[][] Xsi_as,
                                        int maxIt, double tol) {

        // get number of types and time steps
        int ntypes = d.length;
        int m = seq.length;

        // initialize matrix
        double[][][] X0 = new double[m][ntypes][ntypes];
        IntStream.range(0, ntypes)
                .parallel()
                .forEach(i -> {
                    for (int w = 0; w < m; w++) {
                        X0[w][i][i] = rho * (1 - cdf[w][i]);
                    }
                });

        // set up iteration
        double err = 1;
        int it = 0;
        double[][][] X = X0;

        // iterate
        while (err > tol && it < maxIt) {
            double[][][] Xi = new double[m][ntypes][ntypes];

            for (int i = 0; i < ntypes; i++) {
                for (int j = 0; j < ntypes; j++) {
                    // get vectors for convolution
                    double[] y = new double[m];
                    for (int w = 0; w < m; w++) { // multiply elementwise on times
                        for (int k = 0; k < ntypes; k++) { // sum over all types k
                            y[w] += Xsi_s[i][k] * P0[w][k] * X[w][k][j] +
                                    0.5 * Xsi_as[i][k] * (P0[w][i] * X[w][k][j] + P0[w][k] * X[w][i][j]);
                        }
                    }

                    // extract column from the pdf matrix
                    Complex[] Ft = new Complex[m * 2];
                    for (int w = 0; w < m * 2; w++) {
                        Ft[w] = pdfFFT[w][i];
                    }

                    // partially convolve
                    double[] I = convolveFFT(Ft, y, m, dx);

                    // sum
                    for (int w = 0; w < m; w++) {
                        Xi[w][i][j] = X0[w][i][j] + 2 * (1 - d[i]) * I[w];
                    }
                }
            }

            // compute error
            err = getMatrixError3D(X, Xi);

            // update
            X = Xi;
            it++;
        }

        if (it == maxIt) {
            System.err.printf("calcP1 Warning: max iterations reached with error: %.2f%n", err);
        }

        return X;
    }
    

    // Function for calculating branch probabilities
    public static double[][][] calcMTB(double[] a, double[] b, double[] d, double[][] Xsi_s, double[][] Xsi_as,
                                       double[] s, double[] e, HashMap<Integer, UnivariateFunction> P0Map,
                                       int maxIt, double tol, int m) {

        // get number of types and branches
        int ntypes = a.length;
        int nbranches = s.length;
        assert e.length == nbranches; // for each branch, a start and end time must be given

        double[][][] B = new double[nbranches][ntypes][ntypes];

        // for each branch:
        IntStream.range(0, nbranches)
                .parallel()
                .forEach(x -> {
                    double sx = s[x]; // start of the branch
                    double ex = e[x]; // end of the branch
                    assert sx < ex;

                    // generate linearly spaced values between start and end
                    double[] seq = new double[m];
                    double[] age_seq = new double[m];
                    double dx = (ex - sx) / m;
                    for (int w = 0; w < m; w++) {
                        seq[w] = sx + dx * (w + 1);
                        age_seq[w] = seq[w] - sx;
                    }
                    assert seq[m - 1] == ex;

                    // calculate the FFT from PDF of the gamma distribution and interpolate P0 per type, initialize matrix
                    Complex[][] pdfFFT = new Complex[m*2][ntypes];
                    double[][] P0 = new double[m][ntypes];
                    double[][][] X0 = new double[m][ntypes][ntypes];

                    for (int i = 0; i < ntypes; i++) {
                        double[] pdf = new double[m];
                        GammaDistribution gammaDist = new GammaDistribution(b[i], a[i]);
                        for (int w = 0; w < m; w++) {
                            pdf[w] = Math.exp(gammaDist.logDensity(age_seq[w])); // get density
                            P0[w][i] = P0Map.get(i).value(seq[w]); // extrapolate P0
                            X0[w][i][i] = (1 - d[i]) * pdf[w]; // initialize matrix
                        }

                        // perform FFT
                        Complex[] Ft = fft.transform(padZeros(pdf), TransformType.FORWARD);
                        for (int w = 0; w < m*2; w++) {
                            pdfFFT[w][i] = Ft[w];
                        }
                    }

                    // set up iteration
                    double err = 1;
                    int it = 0;
                    double[][][] X = X0;

                    // iterate
                    while (err > tol && it < maxIt) {
                        double[][][] Xi = new double[m][ntypes][ntypes];

                        for (int i = 0; i < ntypes; i++) {
                            for (int j = 0; j < ntypes; j++) {
                                // get vectors for convolution
                                double[] y = new double[m];
                                for (int w = 0; w < m; w++) { // multiply elementwise on times
                                    for (int k = 0; k < ntypes; k++) { // sum over all types k
                                        y[w] += Xsi_s[i][k] * P0[w][k] * X[w][k][j] +
                                                0.5 * Xsi_as[i][k] * (P0[w][i] * X[w][k][j] + P0[w][k] * X[w][i][j]);
                                    }
                                }

                                // extract column from the pdf matrix
                                Complex[] Ft = new Complex[m * 2];
                                for (int w = 0; w < m * 2; w++) {
                                    Ft[w] = pdfFFT[w][i];
                                }

                                // partially convolve
                                double[] I = convolveFFT(Ft, y, m, dx);

                                // sum
                                for (int w = 0; w < m; w++) {
                                    Xi[w][i][j] = X0[w][i][j] + 2 * (1 - d[i]) * I[w];
                                }
                            }
                        }

                        // compute error
                        err = getMatrixError3D(X, Xi);

                        // update
                        X = Xi;
                        it++;
                    }

                    if (it == maxIt) {
                        System.err.printf("calcB branch %d Warning: max iterations reached with error: %.2f%n", x, err);
                    }

                    // take time slice only
                    for (int i = 0; i < ntypes; i++) {
                        for (int j = 0; j < ntypes; j++) {
                            B[x][i][j] = X[m - 1][i][j];
                        }
                    }
                });

        return B;
    }


    // Helper method for calculating the 1-(Manhattan)-distance between two matrices,
    // which is the maximum absolute column sum of the matrices substracted element-wise
    private static double getMatrixError(double[][] X, double[][] Y) {
        int n = X.length; // number of rows
        int m = X[0].length;  // number of columns

        double maxColSum = 0;

        // iterate through each column
        for (int col = 0; col < m; col++) {
            double colSum = 0;

            // sum the absolute differences for this column
            for (int row = 0; row < n; row++) {
                colSum += Math.abs(X[row][col] - Y[row][col]);
            }

            // Update maxColumnSum if this column has a larger sum
            if (colSum > maxColSum) {
                maxColSum = colSum;
            }
        }

        return maxColSum;
    }


    // Helper method for calculating the distance between two 3D arrays
    private static double getMatrixError3D(double[][][] X, double[][][] Y){
        int n = X.length; // number of arrays
        int m = X[0].length;  // number of rows
        int o = X[0][0].length;  // number of columns

        double maxSum = 0;

        // iterate over arrays
        for (int i = 0; i < n; i++) {
            double sum = 0;

            // substract matrices element-wise
            for (int row = 0; row < m; row++) {
                for (int col = 0; col < o; col++) {
                    // sum the absolute differences for each element
                    sum += Math.abs(X[i][row][col] - Y[i][row][col]);
                }
            }

            // update maxSum if this array has a larger sum
            if (sum > maxSum) {
                maxSum = sum;
            }
        }

        return maxSum;
    }

}



