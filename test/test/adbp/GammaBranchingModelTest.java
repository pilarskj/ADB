package test.adbp;

import adbp.GammaBranchingModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.util.Arrays;

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
        model.setInputValue("scale", new RealParameter("3")); //1
        model.setInputValue("shape", new RealParameter("1"));
        model.setInputValue("rho", new RealParameter("0.9")); //0.1
        model.setInputValue("origin", new RealParameter("15")); //12

        model.initAndValidate();

        // calculate tree LL and measure time
        long startTime = System.currentTimeMillis();
        double logL = model.calculateTreeLogLikelihood(tree);
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        System.out.println(logL);
        //System.out.println(elapsedTime);
        //assertEquals(logL, -26.83856, 0.00001);

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
        model.setInputValue("shape", new RealParameter("1"));
        model.setInputValue("rho", new RealParameter("0.6"));
        model.setInputValue("origin", new RealParameter("15"));

        // loop over different scales
        double start = 1;
        double end = 10;
        double step = 0.5;
        FileWriter writer = new FileWriter("examples/logL.txt");
        for (double i = start; i <= end; i += step) {
            model.setInputValue("scale", new RealParameter(Double.toString(i)));
            model.initAndValidate();
            double logL = model.calculateTreeLogLikelihood(tree);
            writer.write(Double.toString(logL));
            writer.write("\n");
        }
        writer.close();
    }
}