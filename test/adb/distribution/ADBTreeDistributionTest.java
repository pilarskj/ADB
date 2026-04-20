package adb.distribution;

import adb.archive.MTBranchingModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

public class ADBTreeDistributionTest {

    // test model parameters and settings
    @Test
    public void testSingleType() throws Exception {
        Parameterization model = new Parameterization();
        model.initByName("nTypes", 1,
                "lifetime", new RealParameter("5"),
                "shape", new IntegerParameter("5"),
                "death", new RealParameter("0.1"),
                "sampling", new RealParameter("0.1"),
                "originTime", 15.0);

        Tree tree = new TreeParser("((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;", false);

        ADBTreeDistribution distribution = new ADBTreeDistribution();
        distribution.initByName("tree", tree,
                "parameterization", model,
                "approx", true,
                "conditionOnOrigin", true);

        double logL = distribution.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
    }

    @Test
    public void testMultiType() throws Exception {

        // define tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/mtADB/profiling/tree_1_size.newick", "IsLabelledNewick", true, "adjustTipHeights", true);

        // initialize and set parameters
        Parameterization model = new Parameterization();
        model.initByName("nTypes", 2,
                "lifetime", new RealParameter("2 5"),
                "shape", new IntegerParameter("1 1"),
                "death", new RealParameter("0.1 0.2"),
                "symTransitions", new RealParameter("0.2 0 0 0.6"),
                "asymTransitions", new RealParameter("0 0.8 0.4 0"),
                "sampling", new RealParameter("0.5 0.5"),
                "originTime", 20.0,
                "originType", 0);

        // calculate tree log-likelihood
        ADBTreeDistribution distribution = new ADBTreeDistribution();
        distribution.initByName("tree", tree,
                "parameterization", model,
                "approx", false,
                "conditionOnOrigin", true);

        double logL = distribution.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
    }

}
