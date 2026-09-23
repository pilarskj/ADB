package adb;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.PositiveInt;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.IntScalar;
import beast.base.spec.type.RealScalar;
import feast.fileio.TreeFromNewickFile;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


// Profile likelihood calculation under ADB
public class ADBProfiling {

    void main() throws Exception {
        profileLogLikCalc();
        // calcLogLikGrid();
        // errorLargeDeath();
        // errorSampling();
    }

    // compare likelihood calculation with exact vs. approximated branch probabilities for trees of different sizes (runtime and error)
    // see https://github.com/pilarskj/ADB-analysis/treesize_comparison
    void profileLogLikCalc() throws Exception {

        String inDir = "/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/";
        String outDir = "src/test/out/";

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader(inDir + "tree_data.csv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(outDir + "tree_likelihood.csv"));
        br.readLine(); // skip header in the input
        bw.write("tree,ntips,exact_time,exact_logL,approx_time,approx_logL\n");  // add header in the output

        // define parameters
        RealScalar<PositiveReal> lifetimeParam = new RealScalarParam<>(10.0, PositiveReal.INSTANCE);
        RealScalar<PositiveReal> shapeRealParam = new RealScalarParam<>(5.0, PositiveReal.INSTANCE);
        IntScalar<PositiveInt> shapeIntParam = new IntScalarParam<>(5, PositiveInt.INSTANCE);
        RealScalar<UnitInterval> deathprobParam = new RealScalarParam<>(0.1, UnitInterval.INSTANCE);
        RealScalar<UnitInterval> rhoParam = new RealScalarParam<>(0.1, UnitInterval.INSTANCE);

        // read each line from the csv
        String line;
        while ((line = br.readLine()) != null) {
            System.out.println(line);
            String[] values = line.split(","); // split the line by comma

            // read values from columns and tree from file
            int treeNr = Integer.parseInt(values[0]);
            int nTips = Integer.parseInt(values[1]);
            double origin = Double.parseDouble(values[2]);
            Tree tree = new TreeFromNewickFile();
            tree.initByName("fileName", inDir + "trees/tree_n" + nTips + "_" + treeNr + ".newick",
                    "IsLabelledNewick", true, "adjustTipHeights", true);

            // record calculation times
            long startTime;
            long endTime;

            GammaBranchingModel exact = new GammaBranchingModel();
            exact.initByName("tree", tree, "origin", origin,
                    "lifetime", lifetimeParam, "shapeReal", shapeRealParam, "deathprob", deathprobParam, "rho", rhoParam,
                    "approx", false);
            startTime = System.nanoTime();
            double exactLik = exact.calculateTreeLogLikelihood(tree);
            endTime = System.nanoTime();
            double exactTime = (endTime - startTime) / 1e+6 ;  // divide to ms

            GammaBranchingModel approx = new GammaBranchingModel();
            approx.initByName("tree", tree, "origin", origin,
                    "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                    "approx", true);
            startTime = System.nanoTime();
            double approxLik = approx.calculateTreeLogLikelihood(tree);
            endTime = System.nanoTime();
            double approxTime = (endTime - startTime) / 1e+6 ;  // divide to ms

            bw.write(treeNr + "," + nTips + "," + exactTime + "," + exactLik + "," + approxTime + "," + approxLik);
            bw.newLine();
        }
        bw.close();
    }


    // compare likelihood calculation on large trees with exact vs. approximated branch probabilities for a range of death probabilities
    // see https://github.com/pilarskj/ADB-analysis/tree/main/treesize_comparison
    void errorLargeDeath() throws Exception {

        String inDir = "/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/";
        String outDir = "src/test/out/";

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader(inDir + "tree_data.csv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(outDir + "tree_likelihood_death.csv"));
        br.readLine(); // skip header in the input
        DecimalFormat df = new DecimalFormat("0.000"); // format death probabilities

        // define parameters
        RealScalar<PositiveReal> lifetimeParam = new RealScalarParam<>(10.0, PositiveReal.INSTANCE);
        RealScalar<PositiveReal> shapeRealParam = new RealScalarParam<>(5.0, PositiveReal.INSTANCE);
        IntScalar<PositiveInt> shapeIntParam = new IntScalarParam<>(5, PositiveInt.INSTANCE);
        RealScalar<UnitInterval> rhoParam = new RealScalarParam<>(0.1, UnitInterval.INSTANCE);

        // read and filter lines in memory
        List<String> filteredLines;
        filteredLines = br.lines()
                .map(String::trim)
                .filter(line -> {
                    String[] values = line.split(",");
                    int treeNr = Integer.parseInt(values[0]);
                    int nTips = Integer.parseInt(values[1]);
                    return (treeNr == 2 || treeNr == 3 || treeNr == 4) && nTips == 5000;})
                .toList();

        // process each tree line in parallel
        List<String> results = filteredLines.parallelStream()
                .flatMap(line -> {
                    System.out.println(line);
                    String[] values = line.split(",");
                    int treeNr = Integer.parseInt(values[0]);
                    int nTips = Integer.parseInt(values[1]);
                    double origin = Double.parseDouble(values[2]);

                    Tree tree = new TreeFromNewickFile();
                    tree.initByName("fileName", inDir + "trees/tree_n" + nTips + "_" + treeNr + ".newick",
                            "IsLabelledNewick", true, "adjustTipHeights", true);

                    return java.util.stream.IntStream.rangeClosed(0, 9)  // define steps
                            .mapToObj(k -> {
                                double deathprob = 0.055 + 0.01 * k;
                                RealScalar<UnitInterval> deathprobParam = new RealScalarParam<>(deathprob, UnitInterval.INSTANCE);

                                GammaBranchingModel exact = new GammaBranchingModel();
                                exact.initByName("tree", tree, "origin", origin,
                                        "lifetime", lifetimeParam, "shapeReal", shapeRealParam, "deathprob", deathprobParam, "rho", rhoParam,
                                        "approx", false);
                                double exactLik = exact.calculateTreeLogLikelihood(tree);

                                GammaBranchingModel approx = new GammaBranchingModel();
                                approx.initByName("tree", tree, "origin", origin,
                                        "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                                        "approx", true);
                                double approxLik = approx.calculateTreeLogLikelihood(tree);

                                return treeNr + "," + df.format(deathprob) + "," + exactLik + "," + approxLik;
                            });
                })
                .toList();

        // write all results to file
        bw.write("tree,deathprob,exact_logL,approx_logL\n");
        for (String resultLine : results) {
            bw.write(resultLine);
            bw.newLine();
        }
        bw.close();
    }


    // compare likelihood calculation with exact vs. approximated branch probabilities for a range of parameter values in a more systematic way
    // see https://github.com/pilarskj/ADB-analysis/tree/main/accuracy_evaluation
    void calcLogLikGrid() throws Exception {

        String inDir = "/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/";
        String outDir = "src/test/out/";

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader(inDir + "trees/trees_grid_sampling.tsv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter(outDir + "loglik/loglik_grid_sampling.csv"));
        br.readLine(); // skip header in the input
        DecimalFormat df = new DecimalFormat("0.00"); // format parameter values

        // read and filter lines in memory
        List<String> lines;
        lines = br.lines().map(String::trim).toList();

        // process each tree line in parallel
        List<String> results = lines.parallelStream()
                .flatMap(line -> {
                    System.out.println(line);
                    String[] values = line.split("\t");
                    int treeNr = Integer.parseInt(values[0]);
                    int shape = Integer.parseInt(values[1]);
                    double rho = Double.parseDouble(values[2]);
                    String treeS = values[3];
                    double origin = Double.parseDouble(values[4]);

                    RealScalar<PositiveReal> lifetimeParam = new RealScalarParam<>(10.0, PositiveReal.INSTANCE);
                    IntScalar<PositiveInt> shapeIntParam = new IntScalarParam<>(shape, PositiveInt.INSTANCE);
                    RealScalar<UnitInterval> rhoParam = new RealScalarParam<>(rho, UnitInterval.INSTANCE);

                    Tree tree = new TreeParser(treeS, true);

                    return java.util.stream.IntStream.rangeClosed(0, 50)
                            .mapToObj(i -> {
                                double deathprob = 0.01 * i;
                                RealScalar<UnitInterval> deathprobParam = new RealScalarParam<>(deathprob, UnitInterval.INSTANCE);
                                // double lifetime = 1 + 0.2 * i;

                                // step sizes
                                GammaBranchingModel ss10 = new GammaBranchingModel();
                                ss10.initByName("tree", tree, "origin", origin,
                                        "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                                        "approx", true, "stepSizeP", (int)Math.pow(2, 10));
                                double ss10Lik = ss10.calculateTreeLogLikelihood(tree);

                                GammaBranchingModel ss12 = new GammaBranchingModel();
                                ss12.initByName("tree", tree, "origin", origin,
                                        "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                                        "approx", true, "stepSizeP", (int)Math.pow(2, 12));
                                double ss12Lik = ss12.calculateTreeLogLikelihood(tree);

                                GammaBranchingModel ss14 = new GammaBranchingModel();
                                ss14.initByName("tree", tree, "origin", origin,
                                        "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                                        "approx", true, "stepSizeP", (int)Math.pow(2, 14));
                                double ss14Lik = ss14.calculateTreeLogLikelihood(tree);

                                return treeNr + ",deathprob," + df.format(deathprob) + "," + ss10Lik + "," + ss12Lik + "," + ss14Lik;

                                /* // exact
                                GammaBranchingModel exact = new GammaBranchingModel();
                                exact.initByName("tree", tree, "origin", origin,
                                        "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                                        "approx", false);
                                double exactLik = exact.calculateTreeLogLikelihood(tree);

                                // approx
                                GammaBranchingModel approx = new GammaBranchingModel();
                                approx.initByName("tree", tree, "origin", origin,
                                        "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                                        "approx", true);
                                double approxLik = approx.calculateTreeLogLikelihood(tree);

                                return treeNr + "," + df.format(deathProb) + "," + exactLik + "," + approxLik;
                                 */
                            });
                })
                .toList();

        // write all results to file
        // bw.write("tree,deathprob,exact_logL,approx_logL\n");
        bw.write("tree,param,value,ss10_logL,ss12_logL,ss14_logL\n");
        for (String resultLine : results) {
            bw.write(resultLine);
            bw.newLine();
        }
        bw.close();
    }


    // quantify error in log-likelihood curves wrt. sampling proportion
    // see https://github.com/pilarskj/ADB-analysis/tree/main/accuracy_evaluation
    void errorSampling() throws Exception {

        String inDir = "/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/";
        String outDir = "src/test/out/";

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader(inDir + "trees/trees_grid_size.tsv")); // or systematic
        BufferedWriter bwLik = new BufferedWriter(new FileWriter(outDir + "/loglik/loglik_grid_size.csv"));
        // BufferedWriter bwTime = new BufferedWriter(new FileWriter(outDir + "/loglik/runtime_systematic_sampling.csv"));
        br.readLine(); // skip header in the input
        String header = "tree,deathprob,stepSize10,stepSize12,stepSize14,stepSize16,exact,BD\n"; // add header in the outputs (or lifetime)
        bwLik.write(header);
        // bwTime.write(header);
        DecimalFormat df = new DecimalFormat("0.00"); // format parameter values

        // range of parameter values
        double min = 0;
        double max = 0.3;
        double step = 0.02;

        // powers for step size
        int[] powers = {10, 12, 14, 16};
        int n = powers.length;

        // read each line from the csv
        String line;
        while ((line = br.readLine()) != null) {
            System.out.println(line);
            String[] values = line.split("\t");
            int treeNr = Integer.parseInt(values[0]);
            int shape = Integer.parseInt(values[1]);
            // double rho = Double.parseDouble(values[2]);
            String treeS = values[3];
            double origin = Double.parseDouble(values[4]);
            double lifetime = 5.0;
            double rho = 0.1;

            Tree tree = new TreeParser(treeS, true);

            // time trackers
            long startTime;
            long endTime;

            // loop over parameter values
            for (double x = min; x <= max; x += step) {
                // double lifetime = x;
                double deathprob = x;

                // parameters
                RealScalar<PositiveReal> lifetimeParam = new RealScalarParam<>(lifetime, PositiveReal.INSTANCE);
                IntScalar<PositiveInt> shapeIntParam = new IntScalarParam<>(shape, PositiveInt.INSTANCE);
                RealScalar<UnitInterval> deathprobParam = new RealScalarParam<>(deathprob, UnitInterval.INSTANCE);
                RealScalar<UnitInterval> rhoParam = new RealScalarParam<>(rho, UnitInterval.INSTANCE);

                // outputs
                double[] loglik = new double[n];
                double[] runtime = new double[n];
                for (int i = 0; i < n; i++) {

                    // ADB calculation with varying step size
                    GammaBranchingModel adb = new GammaBranchingModel();
                    adb.initByName("tree", tree, "origin", origin,
                            "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                            "approx", true, "stepSizeP", (int)Math.pow(2, powers[i]),
                            "useAnalyticalBDSolution", false);
                    startTime = System.nanoTime();
                    loglik[i] = adb.calculateTreeLogLikelihood(tree);
                    endTime = System.nanoTime();
                    runtime[i] = (endTime - startTime) / 1e+6;
                }

                // reference exact calculation
                GammaBranchingModel exact = new GammaBranchingModel();
                exact.initByName("tree", tree, "origin", origin,
                        "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                        "approx", false, "useAnalyticalBDSolution", false);
                double loglikExact = exact.calculateTreeLogLikelihood(tree);

                // reference BD calculation
                double loglikBD;
                if (shape == 1) {
                    GammaBranchingModel bd = new GammaBranchingModel();
                    bd.initByName("tree", tree, "origin", origin,
                            "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                            "useAnalyticalBDSolution", true);
                    startTime = System.nanoTime();
                    loglikBD = bd.calculateTreeLogLikelihood(tree);
                    endTime = System.nanoTime();
                    double runtimeBD = (endTime - startTime) / 1e+6;
                }
                else {loglikBD = Double.NEGATIVE_INFINITY;}

                // write results to files
                String loglikADB = Arrays.stream(loglik).mapToObj(String::valueOf).collect(Collectors.joining(","));
                // String runtimeADB = Arrays.stream(runtime).mapToObj(String::valueOf).collect(Collectors.joining(","));
                String lineLik = treeNr + "," + df.format(x) + "," + loglikADB + "," + loglikExact + "," + loglikBD + "\n";
                // String lineTime = treeNr + "," + df.format(x) + "," + runtimeADB + "," + runtimeBD + "\n";
                bwLik.write(lineLik);
                // bwTime.write(lineTime);
            }
        }
        bwLik.close();
        // bwTime.close();
    }

}