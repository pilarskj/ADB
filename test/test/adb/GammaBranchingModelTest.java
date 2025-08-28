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
import java.util.List;


// Test likelihood calculation under ADB
public class GammaBranchingModelTest {

    // test traversal of tree (for conditioning on the root etc.)
    @Test
    public void testTreeTraversal() throws Exception {

        // define tree
        Tree tree = new TreeParser("((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;", false); // (very small)

        System.out.println("Nodes:");
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            System.out.println(node);
        }

        Node root = tree.getRoot();
        double rootHeight = root.getHeight();
        System.out.println("Root at: " + rootHeight + " (tree height)");

        System.out.println("Left subtree:");
        Node childLeft = root.getChild(0);
        Node childRight = root.getChild(1);
        Node[] nodesLeft = childLeft.getAllChildNodesAndSelf().toArray(new Node[0]); // left subtree
        for (Node node : nodesLeft) {
            System.out.println(node);
        }

        System.out.println("Right subtree:");
        Node[] nodesRight = childRight.getAllChildNodesAndSelf().toArray(new Node[0]); // right subtree
        for (Node node : nodesRight) {
            System.out.println(node);
        }
    }


    // test model parameters and settings
    @Test
    public void testGammaBranchingModel() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // define tree
        Tree tree = new TreeParser("((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;", false); // (very small)
        /* // alternatively, read a tree from .newick file
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "test_data/treeBD.newick", "IsLabelledNewick", true, "adjustTipHeights", true); */
        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("lifetime", new RealParameter("5"));
        model.setInputValue("shapeInteger", new IntegerParameter("5"));
        // model.setInputValue("shapeReal", new RealParameter("5"));
        model.setInputValue("deathprob", new RealParameter("0.1"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("15"));
        // model.setInputValue("maxIterations", 100);
        // model.setInputValue("toleranceP", 1e-12);
        // model.setInputValue("toleranceB", 1e-6);
        // model.setInputValue("stepSizeP", (int)Math.pow(2, 14));
        // model.setInputValue("stepSizeB", (int)Math.pow(2, 14));
        // model.setInputValue("approx", "true");
        // model.setInputValue("conditionOnRoot", "false");
        // model.setInputValue("useAnalyticalBDSolution", "false");
        model.initAndValidate();

        // calculate tree log-likelihood
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
    }


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
    // see https://github.com/pilarskj/ADB-analysis/treesize_comparison
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

    // compare likelihood calculation with exact vs. approximated branch probabilities for a range of death probabilities in a more systematic way
    // see https://github.com/pilarskj/ADB-analysis/accuracy_evaluation
    @Test
    public void testApproximationGrid() throws Exception {

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/trees_grid_sampling.tsv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/loglik_grid_sampling.csv"));
        br.readLine(); // skip header in the input
        DecimalFormat df = new DecimalFormat("0.00"); // format death probabilities

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

}