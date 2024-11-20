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
        model.setInputValue("scale", new RealParameter("3"));
        model.setInputValue("shape", new IntegerParameter("1"));
        model.setInputValue("deathprob", new RealParameter("0"));
        model.setInputValue("rho", new RealParameter("0.9"));
        model.setInputValue("origin", new RealParameter("15")); //12

        model.initAndValidate();

        // calculate tree LL
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
        assertEquals(-328.7123, logL, 0.01);
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
            model.setInputValue("shape", new IntegerParameter(Double.toString(i)));
            model.initAndValidate();
            double logL = model.calculateTreeLogLikelihood(tree);
            writer.write(Double.toString(logL) + "\n");
        }
        writer.close();
    }
}