package adb.distribution;

import adb.archive.BranchList;
import adb.archive.MTBranchingModel;
import adb.tree.AnnotatedTreeParser;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import static adb.archive.MTLogLikColored.calcMTLogLikColored;

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


    @Test
    public void testAnnotated() throws Exception {

        // simulated tree (branch-typed)
        String newick = "((((1[&type=0]:1.378074977,(11[&type=1]:0.9265429956)[&type=1]:0.4515319812)[&type=1]:0.5929298741)[&type=0]:2.088429726,(((3[&type=0]:0.4513672607)[&type=1]:1.262419189,((6[&type=0]:0.05632331623)[&type=1]:0.7255498933)[&type=1]:0.9319132401)[&type=1]:0.2056271944)[&type=0]:2.140020933)[&type=0]:2.073580344)[&type=0]:1.866985079;";
        AnnotatedTreeParser tree = new AnnotatedTreeParser();
        tree.initByName("newick", newick);

        // initialize and set parameters
        Parameterization model = new Parameterization();
        model.initByName("nTypes", 2,
                "lifetime", new RealParameter("2 1"),
                "shape", new IntegerParameter("100 5"),
                "death", new RealParameter("0.1 0.2"),
                "symTransitions", new RealParameter("0.3 0.3 0.1 0.5"),
                "asymTransitions", new RealParameter("0 0.4 0.4 0"),
                "sampling", new RealParameter("0.5 0.5"),
                "originTime", 8.0,
                "originType", 0);

        // calculate tree log-likelihood
        ADBTreeDistribution distribution = new ADBTreeDistribution();
        distribution.initByName("tree", tree,
                "parameterization", model,
                "approx", false,
                "conditionOnOrigin", true);

        double logL = distribution.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
        // currently different from MTLogLikColored, but maybe correct?
    }



    // test model parameters and settings
    @Test
    public void testSingleTypeAnnotated() throws Exception {

        Tree tree = new TreeParser("((1:43.08956348,17:43.08956348):7.807559093,(8:17.26568023,10:17.26568023):33.63144235):8.415882985;", true);
        AnnotatedTreeParser atree = new AnnotatedTreeParser();
        atree.init("((((((1:6.631381528):8.284092381):12.39196391):15.78212566,((((17:2.640722598):8.177954073):8.816746453):11.94528728):11.50885308):7.807559093,((((8:5.916107016):11.34957321,(10:6.247866922):11.0178133):11.55337734):10.62974778):11.44831723):10.68699444):8.415882985;");

        RealParameter lifetime = new RealParameter("10");
        Parameterization model = new Parameterization();
        model.initByName("nTypes", 1,
                "lifetime", lifetime,
                "shape", new IntegerParameter("20"),
                "death", new RealParameter("0.1"),
                "sampling", new RealParameter("0.1"),
                "originTime", 70.0);

        ADBTreeDistribution distribution = new ADBTreeDistribution();
        distribution.initByName("tree", tree,
                "parameterization", model,
                "approx", false,
                "conditionOnOrigin", true);

        double start = 5.0;
        double end = 15.0;
        double step = 0.1;
        int nsteps = (int) ((end - start) / step);
        double[] lifetimes = new double[nsteps];
        for (int i = 0; i < nsteps; i++) {
            lifetimes[i] = start + i * step;
        }

        double[] logL = new double[nsteps];
        double[] alogL = new double[nsteps];
        for (int i = 0; i < nsteps; i++) {
            lifetime.setValue(0, lifetimes[i]);
            logL[i] = distribution.calculateTreeLogLikelihood(tree);
            alogL[i] = distribution.calculateTreeLogLikelihood(atree);
            System.out.println(lifetimes[i] + "," + logL[i] + "," + alogL[i]);
        }
    }

}
