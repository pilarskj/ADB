package test.adb;

import adb.MTBranchingModel;
import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.stream.Collectors;

public class MTProfiling {

    // profile runtime wrt. tree size
    @Test
    public void profileRuntimeSize() throws Exception {

        // define parameters
        // ADB
        int originType = 0;
        RealParameter lifetime = new RealParameter("2 5");
        RealParameter shape = new RealParameter("1 1");
        RealParameter death = new RealParameter("0.1 0.2");
        RealParameter rho = new RealParameter("0.5");
        RealParameter sTransitions = new RealParameter("0.2 0 0 0.6");
        RealParameter asTransitions = new RealParameter("0 0.8 0.4 0");
        int m14 = (int)Math.pow(2, 14);
        int m12 = (int)Math.pow(2, 12);

        // BDMM-Prime
        int ntypes = 2;
        TypeSet typeSet = new TypeSet(ntypes);
        RealParameter startTypePriorProbs = new RealParameter("1 0");
        SkylineVectorParameter birthRate = new SkylineVectorParameter(null,
                new RealParameter("0.09 0.096"), ntypes);
        SkylineVectorParameter deathRate = new SkylineVectorParameter(null,
                new RealParameter("0.05 0.04"), ntypes);
        SkylineMatrixParameter birthRateAmongDemes = new SkylineMatrixParameter(null,
                new RealParameter("0.36 0.064"), ntypes);
        SkylineMatrixParameter migrationRate = new SkylineMatrixParameter(null, new RealParameter("0"), ntypes);
        SkylineVectorParameter samplingRate = new SkylineVectorParameter(null, new RealParameter("0"), ntypes);
        SkylineVectorParameter removalProb = new SkylineVectorParameter(null, new RealParameter("1"), ntypes);


        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/mtADB/profiling/tree_data_size.tsv"));
        BufferedWriter bwTime = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/runtimes_size.csv"));
        BufferedWriter bwLik = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/loglik_size.csv"));
        br.readLine(); // skip header in the input
        bwTime.write("tree,stepSize12,stepSize14,bdmmPrime\n");  // add header in the output
        bwLik.write("tree,stepSize12,stepSize14,bdmmPrime\n");  // add header in the output

        // read each line from the csv
        String line;
        while ((line = br.readLine()) != null) {
            String[] values = line.split("\t"); // split the line by seperator

            // read values from columns
            int treeNr = Integer.parseInt(values[0]);
            double originTime = Double.parseDouble(values[1]);
            System.out.println(treeNr + ", " + originTime);

            String newick = values[3];
            Tree tree = new TreeParser(newick, true);

            // record calculation times
            long startTime;
            long endTime;

            MTBranchingModel stepSize12 = new MTBranchingModel();
            stepSize12.initByName("tree", tree,
                    "lifetime", lifetime, "shape", shape, "death", death, "rho", rho,
                    "sTransitions", sTransitions, "asTransitions", asTransitions,
                    "originType", originType, "originTime", originTime);
            startTime = System.nanoTime();
            double ss12Lik = stepSize12.calculateTreeLogLikelihood(tree);
            endTime = System.nanoTime();
            double ss12Time = (endTime - startTime) / 1e+6;  // divide to ms

            MTBranchingModel stepSize14 = new MTBranchingModel();
            stepSize14.initByName("tree", tree,
                    "lifetime", lifetime, "shape", shape, "death", death, "rho", rho,
                    "sTransitions", sTransitions, "asTransitions", asTransitions,
                    "originType", originType, "originTime", originTime,
                    "stepSizeP", m14, "stepSizeB", m12);
            startTime = System.nanoTime();
            double ss14Lik = stepSize14.calculateTreeLogLikelihood(tree);
            endTime = System.nanoTime();
            double ss14Time = (endTime - startTime) / 1e+6;

            RealParameter processLength = new RealParameter(Double.toString(originTime));
            Parameterization parameterization = new CanonicalParameterization();
            parameterization.initByName("typeSet", typeSet, "processLength", processLength,
                    "birthRate", birthRate, "deathRate", deathRate, "birthRateAmongDemes", birthRateAmongDemes,
                    "migrationRate", migrationRate, "samplingRate", samplingRate, "removalProb", removalProb,
                    "rhoSampling", new TimedParameter(processLength, rho, ntypes));
            BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
            density.initByName("parameterization", parameterization,
                    "startTypePriorProbs", startTypePriorProbs,
                    "tree", tree, "typeLabel", "type", "parallelize", true);
            startTime = System.nanoTime();
            double bPLik = density.calculateLogP();
            endTime = System.nanoTime();
            double bPTime = (endTime - startTime) / 1e+6;

            bwTime.write(treeNr + "," + ss12Time + "," + ss14Time + "," + bPTime);
            bwTime.newLine();
            bwLik.write(treeNr + "," + ss12Lik + "," + ss14Lik + "," + bPLik);
            bwLik.newLine();
        }
        bwTime.close();
        bwLik.close();
    }


    // profile runtime wrt. number of types
    @Test
    public void profileRuntimeTypes() throws Exception {

        int[] ntypes = new int[15];
        for (int i = 0; i < ntypes.length; i++) {
            ntypes[i] = 2 + i;
        }
        System.out.println(Arrays.toString(ntypes));

        // define parameters
        RealParameter rho = new RealParameter("0.5");

        // ADB
        int originType = 0;
        double r = 0.2; // probability of symmetric transition
        int m14 = (int)Math.pow(2, 14);
        int m12 = (int)Math.pow(2, 12);

        HashMap<Integer, RealParameter> lifetime = new HashMap<>();
        HashMap<Integer, RealParameter> shape = new HashMap<>();
        HashMap<Integer, RealParameter> death = new HashMap<>();
        HashMap<Integer, double[][]> sTransMat = new HashMap<>();
        HashMap<Integer, double[][]> asTransMat = new HashMap<>();
        HashMap<Integer, RealParameter> sTransitions = new HashMap<>();
        HashMap<Integer, RealParameter> asTransitions = new HashMap<>();

        for (int n : ntypes) {
            lifetime.put(n, new RealParameter(String.join(" ", Collections.nCopies(n, "1"))));
            shape.put(n, new RealParameter(String.join(" ", Collections.nCopies(n, "1"))));
            death.put(n, new RealParameter(String.join(" ", Collections.nCopies(n, "0.1"))));
            double[][] smat = new double[n][n];
            double[][] asmat = new double[n][n];
            for (int i = 0; i < n; i++) {
                smat[i][i] = r; // set diagonal entries
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        asmat[i][j] = (1-r)/(n-1);
                    }
                }
            }
            sTransMat.put(n, smat);
            asTransMat.put(n, asmat);
            sTransitions.put(n, new RealParameter(convertMatrix(smat, true)));
            asTransitions.put(n, new RealParameter(convertMatrix(asmat,true)));
        }

        // BDMM-Prime
        HashMap<Integer, TypeSet> typeSet = new HashMap<>();
        HashMap<Integer, RealParameter> startTypePriorProbs = new HashMap<>();
        HashMap<Integer, SkylineVectorParameter> birthRate = new HashMap<>();
        HashMap<Integer, SkylineVectorParameter> deathRate = new HashMap<>();
        HashMap<Integer, SkylineMatrixParameter> birthRateAmongDemes = new HashMap<>();
        HashMap<Integer, SkylineMatrixParameter> migrationRate = new HashMap<>();
        HashMap<Integer, SkylineVectorParameter> samplingRate = new HashMap<>();
        HashMap<Integer, SkylineVectorParameter> removalProb = new HashMap<>();

        for (int n : ntypes) {
            typeSet.put(n, new TypeSet(n));
            startTypePriorProbs.put(n, new RealParameter("1 " + String.join(" ", Collections.nCopies(n-1, "0"))));
            birthRate.put(n, new SkylineVectorParameter(null,
                    new RealParameter(String.join(" ", Collections.nCopies(n, "0.18"))), n)); // (1-d)/l * r
            deathRate.put(n, new SkylineVectorParameter(null,
                    new RealParameter(String.join(" ", Collections.nCopies(n, "0.1"))), n)); // d/l
            double[][] rates = asTransMat.get(n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    rates[i][j] *= 0.9; // multiply by (1-d)/l
                }
            }
            birthRateAmongDemes.put(n, new SkylineMatrixParameter(null,
                    new RealParameter(convertMatrix(rates, false)), n));
            migrationRate.put(n, new SkylineMatrixParameter(null, new RealParameter("0"), n));
            samplingRate.put(n, new SkylineVectorParameter(null, new RealParameter("0"), n));
            removalProb.put(n, new SkylineVectorParameter(null, new RealParameter("0"), n));
        }


        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/mtADB/profiling/tree_data_ntypes.tsv"));
        BufferedWriter bwTime = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/runtimes_ntypes.csv"));
        BufferedWriter bwLik = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/loglik_ntypes.csv"));
        br.readLine(); // skip header in the input
        bwTime.write("tree,stepSize12,stepSize14,bdmmPrime\n");  // add header in the output
        bwLik.write("tree,stepSize12,stepSize14,bdmmPrime\n");  // add header in the output

        // read each line from the csv
        String line;
        //int lineNr = 0;
        //int stopLine = 5;
        while ((line = br.readLine()) != null) { // && lineNr < stopLine) {
            //lineNr++;
            String[] values = line.split("\t"); // split the line by seperator

            // read values from columns
            int treeNr = Integer.parseInt(values[0]);
            int n = Integer.parseInt(values[1]);
            double originTime = Double.parseDouble(values[2]);
            System.out.println(treeNr + ", " + n + ", " + originTime);

            String newick = values[3];
            Tree tree = new TreeParser(newick, true);

            // record calculation times
            long startTime;
            long endTime;

            MTBranchingModel stepSize12 = new MTBranchingModel();
            stepSize12.initByName("tree", tree,
                    "lifetime", lifetime.get(n), "shape", shape.get(n), "death", death.get(n), "rho", rho,
                    "sTransitions", sTransitions.get(n), "asTransitions", asTransitions.get(n),
                    "originType", originType, "originTime", originTime);
            startTime = System.nanoTime();
            double ss12Lik = stepSize12.calculateTreeLogLikelihood(tree);
            endTime = System.nanoTime();
            double ss12Time = (endTime - startTime) / 1e+6;  // divide to ms

            MTBranchingModel stepSize14 = new MTBranchingModel();
            stepSize14.initByName("tree", tree,
                    "lifetime", lifetime.get(n), "shape", shape.get(n), "death", death.get(n), "rho", rho,
                    "sTransitions", sTransitions.get(n), "asTransitions", asTransitions.get(n),
                    "originType", originType, "originTime", originTime,
                    "stepSizeP", m14, "stepSizeB", m12);
            startTime = System.nanoTime();
            double ss14Lik = stepSize14.calculateTreeLogLikelihood(tree);
            endTime = System.nanoTime();
            double ss14Time = (endTime - startTime) / 1e+6;

            RealParameter processLength = new RealParameter(Double.toString(originTime));
            Parameterization parameterization = new CanonicalParameterization();
            parameterization.initByName("typeSet", typeSet.get(n), "processLength", processLength,
                    "birthRate", birthRate.get(n), "deathRate", deathRate.get(n), "birthRateAmongDemes", birthRateAmongDemes.get(n),
                    "migrationRate", migrationRate.get(n), "samplingRate", samplingRate.get(n), "removalProb", removalProb.get(n),
                    "rhoSampling", new TimedParameter(processLength, rho, n));
            BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
            density.initByName("parameterization", parameterization,
                    "startTypePriorProbs", startTypePriorProbs.get(n),
                    "tree", tree, "typeLabel", "type", "parallelize", true);
            startTime = System.nanoTime();
            double bPLik = density.calculateLogP();
            endTime = System.nanoTime();
            double bPTime = (endTime - startTime) / 1e+6;

            bwTime.write(treeNr + "," + ss12Time + "," + ss14Time + "," + bPTime);
            bwTime.newLine();
            bwLik.write(treeNr + "," + ss12Lik + "," + ss14Lik + "," + bPLik);
            bwLik.newLine();
        }
        bwTime.close();
        bwLik.close();
    }


    // profile error wrt. sampling proportion
    @Test
    public void profileErrorSampling() throws Exception {

        // define parameters
        // ADB
        int originType = 0;
        RealParameter shape = new RealParameter("1 1");
        RealParameter death = new RealParameter("0.1 0.2");
        RealParameter sTransitions = new RealParameter("0.2 0 0 0.6");
        RealParameter asTransitions = new RealParameter("0 0.8 0.4 0");

        int m8 = (int)Math.pow(2, 8);
        int m10 = (int)Math.pow(2, 10);
        int m12 = (int)Math.pow(2, 12);
        int m14 = (int)Math.pow(2, 14);
        int m16 = (int)Math.pow(2, 16); // too slow

        // BDMM-Prime (equal for all lifetimes)
        int ntypes = 2;
        double[] d = {0.1, 0.2};
        double[][] mat = {{0.2, 0.8}, {0.4, 0.6}}; // transition probability matrix
        TypeSet typeSet = new TypeSet(ntypes);
        RealParameter startTypePriorProbs = new RealParameter("1 0");
        SkylineMatrixParameter migrationRate = new SkylineMatrixParameter(null, new RealParameter("0"), ntypes);
        SkylineVectorParameter samplingRate = new SkylineVectorParameter(null, new RealParameter("0"), ntypes);
        SkylineVectorParameter removalProb = new SkylineVectorParameter(null, new RealParameter("1"), ntypes);

        // compute parameters for sets of lifetimes
        HashMap<double[], RealParameter> lifetime = new HashMap<>();
        HashMap<double[], SkylineVectorParameter> birthRate = new HashMap<>();
        HashMap<double[], SkylineVectorParameter> deathRate = new HashMap<>();
        HashMap<double[], SkylineMatrixParameter> birthRateAmongDemes = new HashMap<>();
        double min = 0.5;
        double max0 = 5.0;
        double max1 = 10.0;
        double step = 0.5;
        for (double x = min; x <= max0; x += step) {
            for (double y = min; y <= max1; y += step) {
                double[] l = {x, y};
                lifetime.put(l, new RealParameter(convertVector(l)));

                double[] bR = new double[ntypes];
                double[][] bRaD = new double[ntypes][ntypes];
                double[] dR = new double[ntypes];
                for (int i = 0; i < ntypes; i++) {
                    bR[i] = mat[i][i] * (1 - d[i]) / l[i];
                    dR[i] = d[i] / l[i];
                    for (int j = 0; j < ntypes; j++) {
                        if (i != j) { bRaD[i][j] = mat[i][j] * (1 - d[i]) / l[i]; }
                    }
                }
                birthRate.put(l, new SkylineVectorParameter(null,
                        new RealParameter(convertVector(bR)), ntypes));
                deathRate.put(l, new SkylineVectorParameter(null,
                        new RealParameter(convertVector(dR)), ntypes));
                birthRateAmongDemes.put(l, new SkylineMatrixParameter(null,
                        new RealParameter(convertMatrix(bRaD, false)), ntypes));
            }
        }
        System.out.println(lifetime.size()); // number of parameter combinations

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/mtADB/profiling/tree_data_sampling.tsv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/loglik_sampling_lifetimes.csv"));
        br.readLine(); // skip header in the input
        bw.write("tree,lifetime0,lifetime1,stepSize10,stepSize12,stepSize14,bdmmPrime\n");  // add header in the output
        DecimalFormat df = new DecimalFormat("0.0"); // for parameter values

        // read each line from the csv
        String line;
        while ((line = br.readLine()) != null) {
            String[] values = line.split("\t"); // split the line by seperator

            // read values from columns
            int treeNr = Integer.parseInt(values[0]);
            double samplingProb = Double.parseDouble(values[1]);
            RealParameter rho = new RealParameter(Double.toString(samplingProb));
            double originTime = Double.parseDouble(values[2]);
            System.out.println(treeNr + ", " + samplingProb + ", " + originTime);

            String newick = values[3];
            Tree tree = new TreeParser(newick, true);

            for (double[] l : lifetime.keySet()) {
                MTBranchingModel stepSize10 = new MTBranchingModel();
                stepSize10.initByName("tree", tree,
                        "lifetime", lifetime.get(l), "shape", shape, "death", death, "rho", rho,
                        "sTransitions", sTransitions, "asTransitions", asTransitions,
                        "originType", originType, "originTime", originTime,
                        "stepSizeP", m10, "stepSizeB", m8);
                double ss10Lik = stepSize10.calculateTreeLogLikelihood(tree);

                MTBranchingModel stepSize12 = new MTBranchingModel();
                stepSize12.initByName("tree", tree,
                        "lifetime", lifetime.get(l), "shape", shape, "death", death, "rho", rho,
                        "sTransitions", sTransitions, "asTransitions", asTransitions,
                        "originType", originType, "originTime", originTime);
                double ss12Lik = stepSize12.calculateTreeLogLikelihood(tree);

                MTBranchingModel stepSize14 = new MTBranchingModel();
                stepSize14.initByName("tree", tree,
                        "lifetime", lifetime.get(l), "shape", shape, "death", death, "rho", rho,
                        "sTransitions", sTransitions, "asTransitions", asTransitions,
                        "originType", originType, "originTime", originTime,
                        "stepSizeP", m14, "stepSizeB", m12);
                double ss14Lik = stepSize14.calculateTreeLogLikelihood(tree);

                RealParameter processLength = new RealParameter(Double.toString(originTime));
                Parameterization parameterization = new CanonicalParameterization();
                parameterization.initByName("typeSet", typeSet, "processLength", processLength,
                        "birthRate", birthRate.get(l), "deathRate", deathRate.get(l), "birthRateAmongDemes", birthRateAmongDemes.get(l),
                        "migrationRate", migrationRate, "samplingRate", samplingRate, "removalProb", removalProb,
                        "rhoSampling", new TimedParameter(processLength, rho, ntypes));
                BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
                density.initByName("parameterization", parameterization,
                        "startTypePriorProbs", startTypePriorProbs,
                        "tree", tree, "typeLabel", "type", "parallelize", true);
                double bPLik = density.calculateLogP();

                bw.write(treeNr + "," + df.format(l[0]) + "," + df.format(l[1]) + "," + ss10Lik + "," + ss12Lik + "," + ss14Lik + "," + bPLik);
                bw.newLine();
            }
        }
        bw.close();
    }


    // Helper function for converting matrix to String
    private static String convertMatrix(double[][] x, boolean zero) {
        String s;
        if (zero) {
            s = Arrays.stream(x)
                    .flatMapToDouble(Arrays::stream)
                    .mapToObj(Double::toString)
                    .collect(Collectors.joining(" "));
        } else {
            s = Arrays.stream(x)
                    .flatMapToDouble(Arrays::stream)
                    .filter(v -> v != 0.0)
                    .mapToObj(Double::toString)
                    .collect(Collectors.joining(" "));
        }
        return(s);
    }

    // Helper function for converting vector to String
    private static String convertVector(double[] x) {
        String s = Arrays.stream(x)
                .mapToObj(Double::toString)
                .collect(Collectors.joining(" "));
        return(s);
    }

}
