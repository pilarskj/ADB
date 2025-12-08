package test.adb;

import adb.MTBranchingModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

public class MTBranchingModelTest {

    // test model parameters and settings
    @Test
    public void testMTBranchingModel() throws Exception {

        // define tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/mtADB/profiling/tree_1_size.newick", "IsLabelledNewick", true, "adjustTipHeights", true);

        // initialize and set parameters
        MTBranchingModel model = new MTBranchingModel();
        model.initByName("tree", tree,
                "lifetime", new RealParameter("2 5"),
                "shape", new RealParameter("1 1"),
                "death", new RealParameter("0.1 0.2"),
                "sTransitions", new RealParameter("0.2 0 0 0.6"),
                "asTransitions", new RealParameter("0 0.8 0.4 0"),
                "rho", new RealParameter("0.5"),
                "originTime", 20.0,
                "originType", 0);

        // calculate tree log-likelihood
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
    }

}
