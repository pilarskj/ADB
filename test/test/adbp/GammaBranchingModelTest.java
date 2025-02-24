package test.adbp;

import adbp.GammaBranchingModel;
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

import static org.junit.jupiter.api.Assertions.assertEquals;

public class GammaBranchingModelTest {

    @Test
    public void testGammaBranchingModel() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();
        /*
        // define tree
        Tree tree = new TreeParser("((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;", false); // (very small)
        Tree tree = new TreeParser( // tree with very short or long branches
                "(((((0:50,((1:1e-07,2:1e-07)9:1e-07,3:2e-07)10:49.9999998)11:1.000000012e-07,4:50.0000001)12:1.000000012e-07,5:50.0000002)13:1.000000012e-07,(6:1e-07,7:1e-07)14:50.0000002)15:1.000000012e-07,8:50.0000004)16:0;",
                true);
        */
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/test2/rhoLow/trees/tree_1.newick", //"/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/calculations/tree_bd.newick"
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("lifetime", new RealParameter("5"));
        model.setInputValue("shapeInteger", new IntegerParameter("5"));
        //model.setInputValue("shapeReal", new RealParameter("0.1"));
        model.setInputValue("deathprob", new RealParameter("0.05")); //0.1
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("50"));
        model.setInputValue("approx", "true");

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
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/simulation/tree_1.upgma_mod.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);
        model.setInputValue("approx", "true");

        // set parameters
        model.setInputValue("lifetime", new RealParameter("10"));
        model.setInputValue("shapeInteger", new IntegerParameter("10"));
        model.setInputValue("deathprob", new RealParameter("0.5"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("60.5037032882197"));

        //model.setInputValue("toleranceP", 1e-6);
        //model.setInputValue("stepSizeP", (int)Math.pow(2, 18));
        //model.setInputValue("maxIterations", 200);

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
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/simulation/tree_1.upgma_2short.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);
        model.setInputValue("tree", tree);

        // set parameters
        //model.setInputValue("shapeInteger", new IntegerParameter("1"));
        model.setInputValue("lifetime", new RealParameter("5"));
        model.setInputValue("deathprob", new RealParameter("0.1"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("60.5037032882197"));
        model.setInputValue("approx", true);

        //model.setInputValue("stepSizeP", (int)Math.pow(2, 18));

        // loop over different parameter values
        double start = 1;
        double end = 50;
        double step = 1;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/simulation/loglik_shape.csv", true);
        DecimalFormat df = new DecimalFormat("0");
        for (double i = start; i <= end; i += step) {
            model.setInputValue("shapeInteger", new IntegerParameter(Double.toString(i)));
            //model.setInputValue("shapeReal", new RealParameter(Double.toString(i)));
            model.initAndValidate();
            double logL = model.calculateTreeLogLikelihood(tree);
            writer.write( df.format(i) + "," + logL + ",2short\n"); //,adbp_approx
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