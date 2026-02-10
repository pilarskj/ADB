package test.adb;

import adb.GammaBranchingModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


// Test likelihood calculation under ADB
public class ADBProfiling {

    // compare likelihood calculation with exact vs. approximated branch probabilities for trees of different sizes (runtime and error)
    // see https://github.com/pilarskj/ADB-analysis/treesize_comparison
    @Test
    public void testApproximation() throws Exception {

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/tree_data.csv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/tree_likelihood.csv"));
        br.readLine(); // skip header in the input
        bw.write("tree,ntips,exact_time,exact_logL,approx_time,approx_logL\n");  // add header in the output

        // read each line from the csv
        String line;
        while ((line = br.readLine()) != null) {
            System.out.println(line);
            String[] values = line.split(","); // split the line by comma

            // read values from columns and tree from file
            int treeNr = Integer.parseInt(values[0]);
            int nTips = Integer.parseInt(values[1]);
            Tree tree = new TreeFromNewickFile();
            tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/trees/tree_n" + nTips + "_" + treeNr + ".newick",
                    "IsLabelledNewick", true, "adjustTipHeights", true);

            // record calculation times
            long startTime;
            long endTime;

            GammaBranchingModel exact = new GammaBranchingModel();
            exact.initByName("tree", tree,
                    "lifetime", new RealParameter("10"),
                    "shapeInteger", new IntegerParameter("5"),
                    "deathprob", new RealParameter("0.1"),
                    "rho", new RealParameter("0.1"),
                    "origin", new RealParameter(values[2]),
                    "approx", false);
            startTime = System.nanoTime();
            double exactLik = exact.calculateTreeLogLikelihood(tree);
            endTime = System.nanoTime();
            double exactTime = (endTime - startTime) / 1e+6 ;  // divide to ms

            GammaBranchingModel approx = new GammaBranchingModel();
            approx.initByName("tree", tree,
                    "lifetime", new RealParameter("10"),
                    "shapeInteger", new IntegerParameter("5"),
                    "deathprob", new RealParameter("0.1"),
                    "rho", new RealParameter("0.1"),
                    "origin", new RealParameter(values[2]),
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
    @Test
    public void testApproximationSeqParallel() throws Exception {

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/tree_data.csv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/tree_likelihood_death.csv"));
        br.readLine(); // skip header in the input
        DecimalFormat df = new DecimalFormat("0.000"); // format death probabilities

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
                    String originS = values[2];

                    Tree tree = new TreeFromNewickFile();
                    tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/trees/tree_n" + nTips + "_" + treeNr + ".newick",
                            "IsLabelledNewick", true, "adjustTipHeights", true);

                    return java.util.stream.IntStream.rangeClosed(0, 9)  // define steps
                            .mapToObj(k -> {
                                double deathProb = 0.055 + 0.01 * k;
                                // exact
                                GammaBranchingModel exact = new GammaBranchingModel();
                                exact.initByName("tree", tree,
                                        "lifetime", new RealParameter("10"),
                                        "shapeInteger", new IntegerParameter("5"),
                                        "deathprob", new RealParameter(Double.toString(deathProb)),
                                        "rho", new RealParameter("0.1"),
                                        "origin", new RealParameter(originS),
                                        "approx", false);
                                double exactLik = exact.calculateTreeLogLikelihood(tree);

                                // approx
                                GammaBranchingModel approx = new GammaBranchingModel();
                                approx.initByName("tree", tree,
                                        "lifetime", new RealParameter("10"),
                                        "shapeInteger", new IntegerParameter("5"),
                                        "deathprob", new RealParameter(Double.toString(deathProb)),
                                        "rho", new RealParameter("0.1"),
                                        "origin", new RealParameter(originS),
                                        "approx", true);
                                double approxLik = approx.calculateTreeLogLikelihood(tree);

                                return treeNr + "," + df.format(deathProb) + "," + exactLik + "," + approxLik;
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
    @Test
    public void testApproximationGrid() throws Exception {

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/trees/trees_grid_sampling.tsv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/loglik/loglik_grid_sampling.csv"));
        br.readLine(); // skip header in the input
        DecimalFormat df = new DecimalFormat("0.0"); // format parameter values

        // read and filter lines in memory
        List<String> lines;
        lines = br.lines().map(String::trim).toList();

        // process each tree line in parallel
        List<String> results = lines.parallelStream()
                .flatMap(line -> {
                    System.out.println(line);
                    String[] values = line.split("\t");
                    int treeNr = Integer.parseInt(values[0]);
                    String kS = values[1];
                    String rhoS = values[2];
                    String treeS = values[3];
                    String originS = values[4];

                    Tree tree = new TreeParser(treeS, true);

                    return java.util.stream.IntStream.rangeClosed(0, 50)
                            .mapToObj(i -> {
                                double deathProb = 0.01 * i;
                                // double lifetime = 1 + 0.2 * i;

                                // step sizes
                                GammaBranchingModel ss10 = new GammaBranchingModel();
                                ss10.initByName("tree", tree,
                                        "lifetime", new RealParameter("5"), // Double.toString(lifetime)
                                        "shapeInteger", new IntegerParameter(kS),
                                        "deathprob", new RealParameter(Double.toString(deathProb)),
                                        "rho", new RealParameter(rhoS),
                                        "origin", new RealParameter(originS),
                                        "approx", true,
                                        "useAnalyticalBDSolution", false,
                                        "stepSizeP", (int)Math.pow(2, 10)
                                        );
                                double ss10Lik = ss10.calculateTreeLogLikelihood(tree);

                                GammaBranchingModel ss12 = new GammaBranchingModel();
                                ss12.initByName("tree", tree,
                                        "lifetime", new RealParameter("5"), // Double.toString(lifetime)
                                        "shapeInteger", new IntegerParameter(kS),
                                        "deathprob", new RealParameter(Double.toString(deathProb)),
                                        "rho", new RealParameter(rhoS),
                                        "origin", new RealParameter(originS),
                                        "approx", true,
                                        "useAnalyticalBDSolution", false,
                                        "stepSizeP", (int)Math.pow(2, 12)
                                );
                                double ss12Lik = ss12.calculateTreeLogLikelihood(tree);

                                GammaBranchingModel ss14 = new GammaBranchingModel();
                                ss14.initByName("tree", tree,
                                        "lifetime", new RealParameter("5"), // Double.toString(lifetime)
                                        "shapeInteger", new IntegerParameter(kS),
                                        "deathprob", new RealParameter(Double.toString(deathProb)),
                                        "rho", new RealParameter(rhoS),
                                        "origin", new RealParameter(originS),
                                        "approx", true,
                                        "useAnalyticalBDSolution", false,
                                        "stepSizeP", (int)Math.pow(2, 14)
                                );
                                double ss14Lik = ss14.calculateTreeLogLikelihood(tree);

                                /*
                                // exact
                                GammaBranchingModel exact = new GammaBranchingModel();
                                exact.initByName("tree", tree,
                                        "lifetime", new RealParameter("10"),
                                        "shapeInteger", new IntegerParameter(kS),
                                        "deathprob", new RealParameter(Double.toString(deathProb)),
                                        "rho", new RealParameter("0.1"),
                                        "origin", new RealParameter(originS),
                                        "approx", false,
                                        "useAnalyticalBDSolution", false);
                                double exactLik = exact.calculateTreeLogLikelihood(tree);

                                // approx
                                GammaBranchingModel approx = new GammaBranchingModel();
                                approx.initByName("tree", tree,
                                        "lifetime", new RealParameter("10"),
                                        "shapeInteger", new IntegerParameter(kS),
                                        "deathprob", new RealParameter(Double.toString(deathProb)),
                                        "rho", new RealParameter("0.1"),
                                        "origin", new RealParameter(originS),
                                        "approx", true,
                                        "useAnalyticalBDSolution", false);
                                double approxLik = approx.calculateTreeLogLikelihood(tree);
                                */
                                //return treeNr + "," + df.format(deathProb) + "," + exactLik + "," + approxLik;
                                return treeNr + ",deathprob," + df.format(deathProb) + "," + ss10Lik + "," + ss12Lik + "," + ss14Lik;
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
    @Test
    public void profileErrorSampling() throws Exception {

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/trees/trees_grid_size.tsv")); // or systematic
        BufferedWriter bwLik = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/loglik/loglik_grid_size.csv"));
        // BufferedWriter bwTime = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/loglik/runtime_systematic_sampling.csv"));
        br.readLine(); // skip header in the input
        // String header = "tree,lifetime,stepSize10,stepSize12,stepSize14,stepSize16,BD\n"; // add header in the outputs
        String header = "tree,deathprob,stepSize10,stepSize12,stepSize14,stepSize16,exact,BD\n";
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
            String shapeS = values[1];
            // String rhoS = values[2];
            String treeS = values[3];
            String originS = values[4];

            Tree tree = new TreeParser(treeS, true);

            // time trackers
            long startTime;
            long endTime;

            // loop over parameter values
            for (double x = min; x <= max; x += step) {
                // String lifetimeS = Double.toString(x);
                String deathprobS = Double.toString(x);
                // outputs
                double[] runtime = new double[n];
                double[] loglik = new double[n];
                for (int i = 0; i < n; i++) {
                    // ADB calculation with varying step size
                    GammaBranchingModel adb = new GammaBranchingModel();
                    adb.initByName("tree", tree,
                            "lifetime", new RealParameter("5"), // lifetimeS
                            "shapeInteger", new IntegerParameter(shapeS), // shapeS
                            "deathprob", new RealParameter(deathprobS), // "0.1"
                            "rho", new RealParameter("0.1"), // rhoS
                            "origin", new RealParameter(originS),
                            "approx", true,
                            "useAnalyticalBDSolution", false,
                            "stepSizeP", (int) Math.pow(2, powers[i]));
                    startTime = System.nanoTime();
                    loglik[i] = adb.calculateTreeLogLikelihood(tree);
                    endTime = System.nanoTime();
                    runtime[i] = (endTime - startTime) / 1e+6;
                }

                // reference exact calculation
                GammaBranchingModel exact = new GammaBranchingModel();
                exact.initByName("tree", tree,
                        "lifetime", new RealParameter("5"),
                        "shapeInteger", new IntegerParameter(shapeS),
                        "deathprob", new RealParameter(deathprobS),
                        "rho", new RealParameter("0.1"),
                        "origin", new RealParameter(originS),
                        "approx", false,
                        "useAnalyticalBDSolution", false);
                double loglikExact = exact.calculateTreeLogLikelihood(tree);

                // reference BD calculation
                double loglikBD;
                if (shapeS.equals("1")) {
                    GammaBranchingModel bd = new GammaBranchingModel();
                    bd.initByName("tree", tree,
                            "lifetime", new RealParameter("5"),
                            "shapeInteger", new IntegerParameter("1"),
                            "deathprob", new RealParameter(deathprobS),
                            "rho", new RealParameter("0.1"),
                            "origin", new RealParameter(originS),
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
                String lineLik = treeNr + "," + df.format(x) + "," + loglikADB + "," + loglikExact + "," + loglikBD + "\n"; //
                // String lineTime = treeNr + "," + df.format(x) + "," + runtimeADB + "," + runtimeBD + "\n";
                bwLik.write(lineLik);
                // bwTime.write(lineTime);
            }
        }
        bwLik.close();
        // bwTime.close();
    }

}