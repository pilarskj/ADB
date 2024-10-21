package test.adbp;

import adbp.GammaBranchingModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class GammaBranchingModelTest {

    @Test
    public void testGammaBranching() throws Exception {

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
        model.setInputValue("scale", new RealParameter("3"));
        model.setInputValue("shape", new RealParameter("1"));
        model.setInputValue("deathprob", new RealParameter("0"));
        model.setInputValue("rho", new RealParameter("0.9"));
        model.setInputValue("origin", new RealParameter("15")); //12

        model.initAndValidate();

        // calculate tree LL
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
        assertEquals(logL, -328.7123, 0.01);

        // using Streams for calcB (parallel computation per branch), the runtime reduces from >3s to ~950ms!
        // but still, a MCMC chain of 100M samples would run for 3 years...
    }

    @Test
    public void testGammaBranchingSeq() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("scale", new RealParameter("3"));
        model.setInputValue("deathprob", new RealParameter("0"));
        model.setInputValue("rho", new RealParameter("0.9"));
        model.setInputValue("origin", new RealParameter("15"));

        // loop over different scales
        double start = 0.1;
        double end = 5;
        double step = 0.1;
        FileWriter writer = new FileWriter("examples/logL.txt");
        for (double i = start; i <= end; i += step) {
            model.setInputValue("shape", new RealParameter(Double.toString(i)));
            model.initAndValidate();
            double logL = model.calculateTreeLogLikelihood(tree);
            writer.write(Double.toString(logL));
            writer.write("\n");
        }
        writer.close();
    }


    /* for optimizing m (the array size for FFT):
    set mP and mB as parameters of the GammaBranchingModel!
    @Test
    public void testM() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
        //tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/efficiency/tree_b20.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("scale", new RealParameter("3")); //3, 0.1
        model.setInputValue("shape", new RealParameter("0.1")); //1, 0.1, 20
        model.setInputValue("deathprob", new RealParameter("0"));
        model.setInputValue("rho", new RealParameter("0.9"));
        model.setInputValue("origin", new RealParameter("15"));

        // exact LL value for shape = 1 (BD model)
        double LL_exact = -328.7123;
        // LL value for shape = 20 and m = 2^17
        // double LL_exact = -9429.179397226246;
        // LL value for shape = 0.1 and m = 2^17
        // double LL_exact = -845.9254038768404;
        // LL value for shape = 20 and m = 2^17 and tree with true shape = 20
        // double LL_exact = -118.35404596644366;

        // loop over different m
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/efficiency/LL_err_b0.1_mB");
        writer.write("power,logL,err\n");

        for (int i = 11; i <= 17; i += 1) {
            int mP = (int)Math.pow(2, 14);
            int mB = (int)Math.pow(2, i);
            model.setInputValue("mP", new IntegerParameter(Integer.toString(mP)));
            model.setInputValue("mB", new IntegerParameter(Integer.toString(mB)));
            model.initAndValidate();
            // calculate LL
            double logL = model.calculateTreeLogLikelihood(tree);

            // calculate error
            double err = Math.abs(logL - LL_exact);

            // System.out.println(logL);
            writer.write(Integer.toString(i));
            writer.write(",");
            writer.write(Double.toString(logL));
            writer.write(",");
            writer.write(Double.toString(err));
            writer.write("\n");
        }
        writer.close();
    }
     */
}