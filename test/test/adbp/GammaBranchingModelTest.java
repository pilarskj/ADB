package test.adbp;

import adbp.GammaBranchingModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class GammaBranchingModelTest {

    @Test
    public void testGammaBranchingModel() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // define tree (very small)
        // Tree tree = new TreeParser("((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;", false);

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("lifetime", new RealParameter("2")); //3
        model.setInputValue("shape", new IntegerParameter("1"));
        model.setInputValue("deathprob", new RealParameter("0.2")); //0
        model.setInputValue("rho", new RealParameter("0.9"));
        model.setInputValue("origin", new RealParameter("15")); //12
        model.setInputValue("approx", "false");

        model.initAndValidate();

        // calculate tree LL
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
        // assertEquals(-328.7123, logL, 0.01); (before adding the tree factor)
    }


    @Test
    public void testGammaBranchingInference() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/test2/b1/trees/tree_1.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);
        model.setInputValue("approx", "false");

        // set parameters
        model.setInputValue("lifetime", new RealParameter("5"));
        model.setInputValue("shape", new IntegerParameter("1"));
        model.setInputValue("deathprob", new RealParameter("0.1"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("42.8106111505566"));

        model.initAndValidate();

        // calculate tree LL
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
    }


    @Test
    public void testGammaBranchingGrid() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/identifiability/trees/tree_1.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);
        model.setInputValue("approx", "false");
        model.setInputValue("origin", new RealParameter("27.1906287767395"));
        model.setInputValue("rho", new RealParameter("0.1")); // set to true value

        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/identifiability/param_grid.csv"));
        BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/identifiability/param_grid_out.csv"));

        String line;
        boolean isHeader = true;

        // read each line from the csv
        while ((line = br.readLine()) != null) {
            System.out.println(line);
            // split the line by comma
            String[] values = line.split(",");

            // if it's the first line (header), add a new column name
            if (isHeader) {
                isHeader = false;
                bw.write(line + ",logL");  // Add a header for the new column
                bw.newLine();
                continue;
            }

            // read values from columns a, b, and d (assuming they are in known positions)
            // set parameters
            model.setInputValue("scale", new RealParameter(values[0]));
            model.setInputValue("shape", new IntegerParameter(values[1]));
            model.setInputValue("deathprob", new RealParameter(values[2]));
            model.initAndValidate();

            // calculate tree LL
            double logL = model.calculateTreeLogLikelihood(tree);
            bw.write(line + "," + logL);
            bw.newLine();
        }
        bw.close();
    }


    @Test
    public void testGammaBranchingSeq() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/test1/trees/tree_1.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);
        model.setInputValue("tree", tree);

        // set parameters
        // model.setInputValue("shape", new IntegerParameter("1"));
        model.setInputValue("lifetime", new RealParameter("5"));
        model.setInputValue("deathprob", new RealParameter("0.1"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("60.5037032882197")); //17.47451 32.10385 33.59449
        model.setInputValue("approx", true);

        // loop over different parameter values
        double start = 1;
        double end = 50;
        double step = 1;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/test1/tree_1_loglik_before.csv", true);
        DecimalFormat df = new DecimalFormat("0.0");
        for (double i = start; i <= end; i += step) {
            model.setInputValue("shape", new IntegerParameter(Double.toString(i)));
            model.initAndValidate();
            double logL = model.calculateTreeLogLikelihood(tree);
            // add n, param, value, loglik
            writer.write(i + "," + logL + "\n"); //df.format(i)
        }
        writer.close();
    }


    @Test
    public void testApproximationN() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();
        // set parameters
        model.setInputValue("shape", new IntegerParameter("5"));
        model.setInputValue("scale", new RealParameter("0.4"));
        model.setInputValue("deathprob", new RealParameter("0.2"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("approx", true);

        // specify numbers of taxa and tree origins
        int[] ntaxa = {5, 10, 50, 100, 500, 1000, 2500, 5000};
        double[] origin = {16.61848, 19.69262, 27.32094, 26.82863, 46.80538, 39.73703, 44.47180, 41.78479};

        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/approx/approx_ntaxa.csv", true);
        for (int i = 0; i < ntaxa.length; i++)  {
            int n = ntaxa[i];
            double o = origin[i];
            // get tree
            Tree tree = new TreeFromNewickFile();
            tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/approx/tree_" + n + "_tips.newick",
                    "IsLabelledNewick", true,
                    "adjustTipHeights", true);
            model.setInputValue("tree", tree);
            model.setInputValue("origin", new RealParameter(Double.toString(o)));
            model.initAndValidate();

            // record execution time and log likelihood
            long startTime = System.nanoTime();
            double logL = model.calculateTreeLogLikelihood(tree);
            long endTime = System.nanoTime();
            double duration = (endTime - startTime) / 1e+6;  // divide to ms
            writer.write(n + "," + logL + "," + duration + ",true" +"\n");
        }
        writer.close();
    }


    @Test
    public void testApproximationO() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();
        // set parameters
        model.setInputValue("shape", new IntegerParameter("4"));
        model.setInputValue("deathprob", new RealParameter("0.2"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("approx", true);

        // specify numbers of taxa and tree origins
        double[] a_seq = {0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5};
        double[] origin = {26.3021108427124, 54.4472552478273, 82.2056146691726, 90.9562602553851, 124.761235895773, 144.422837524865, 180.716604472744, 247.852045064087, 315.627710546529, 269.205117746096};

        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/approx/approx_scale.csv", true);
        for (int i = 0; i < a_seq.length; i++)  {
            double a = a_seq[i];
            double o = origin[i];

            // get tree
            Tree tree = new TreeFromNewickFile();
            tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/approx/tree_scale_" + a + ".newick",
                    "IsLabelledNewick", true,
                    "adjustTipHeights", true);
            model.setInputValue("tree", tree);
            model.setInputValue("origin", new RealParameter(Double.toString(o)));
            model.setInputValue("scale", new RealParameter(Double.toString(a)));
            model.initAndValidate();

            // record execution time and log likelihood
            long startTime = System.nanoTime();
            double logL = model.calculateTreeLogLikelihood(tree);
            long endTime = System.nanoTime();
            double duration = (endTime - startTime) / 1e+6;  // divide to ms
            writer.write(a + "," + logL + "," + duration + ",true" +"\n");
        }
        writer.close();
    }

    @Test
    public void testApproximationP() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/approx/tree_scale_0.5.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);
        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("shape", new IntegerParameter("4"));
        model.setInputValue("scale", new RealParameter("0.5"));
        model.setInputValue("deathprob", new RealParameter("0.2"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("26.3021108427124"));
        model.setInputValue("approx", true);

        // loop over different parameter values
        double start = 0.1;
        double end = 1;
        double step = 0.1;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/approx/approx_params.csv", true);
        DecimalFormat df = new DecimalFormat("0.0");
        for (double i = start; i <= end; i += step) {
            model.setInputValue("rho", new RealParameter(Double.toString(i)));
            model.initAndValidate();
            double logL = model.calculateTreeLogLikelihood(tree);
            // add param, value, loglik, approx
            writer.write("rho," + df.format(i) + "," + logL + ",true\n");
        }
        writer.close();
    }
}