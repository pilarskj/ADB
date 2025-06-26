package test.adbp;

import adbp.GammaBranchingModel;
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
        tree.initByName("fileName", "test_data/example.newick", "IsLabelledNewick", true, "adjustTipHeights", true); */
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


    // compare runtime for exact vs. approximated branch probabilities
    // see https://github.com/pilarskj/ADB-analysis/treesize_comparison
    @Test
    public void testApproximationRuntime() throws Exception {

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/tree_data.csv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/treesize_comparison/runtime_likelihood.csv"));
        br.readLine(); // skip header in the input
        bw.write("tree,ntips,exact_time,approx_time\n");  // add header in the output

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
            exact.calculateTreeLogLikelihood(tree);
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
            approx.calculateTreeLogLikelihood(tree);
            endTime = System.nanoTime();
            double approxTime = (endTime - startTime) / 1e+6 ;  // divide to ms

            bw.write(treeNr + "," + nTips + "," + exactTime + "," + approxTime);
            bw.newLine();
        }
        bw.close();
    }

}