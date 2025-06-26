package test.adbp;

import adbp.GammaBranchingModel;
import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import bdsky.evolution.speciation.BirthDeathSkylineModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.text.DecimalFormat;

import static org.junit.jupiter.api.Assertions.assertEquals;


// Comparison of likelihood calculation in ADB, BDSKY and BDMM-Prime for tree with shape = 1 (BD case)
public class BDTest {

    // compare on an example tree
    @Test
    public void testBDSKY() throws Exception {
        // initialize
        BirthDeathSkylineModel model = new BirthDeathSkylineModel();

        // get tree with 100 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "test_data/treeBD.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("birthRate", new RealParameter("0.18"));
        model.setInputValue("deathRate", new RealParameter("0.02"));
        model.setInputValue("samplingRate", new RealParameter("0"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("50"));
        model.setInputValue("contemp", "true");
        model.setInputValue("conditionOnSurvival", "true");
        //model.setInputValue("conditionOnRoot", "true");

        model.initAndValidate();

        // calculate tree LL
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL); // -348.3790 value matches with ADB without tree factor (has been discarded)
    }

    @Test
    public void testBDMM() throws Exception {

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "test_data/treeBD.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // initialize
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new RealParameter("50"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.18")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.02")),
                "birthRateAmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("50"),
                        new RealParameter("0.1")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", new RealParameter("1"),
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
                // "conditionOnRoot", true,
                "useAnalyticalSingleTypeSolution", true);

        double logL = density.calculateLogP();
        System.out.println(logL); // -643.4968 value matches with ADB
    }

    @Test
    public void testADB() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "test_data/treeBD.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);
        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("lifetime", new RealParameter("5"));
        model.setInputValue("shapeInteger", new IntegerParameter("1"));
        model.setInputValue("deathprob", new RealParameter("0.1"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("50"));
        model.setInputValue("approx", "approx");

        model.initAndValidate();

        // calculate tree LL
        model.setInputValue("useAnalyticalBDSolution", "true");
        double bd = model.calculateTreeLogLikelihood(tree);
        System.out.println("BD logL = " + bd); // -643.4968

        model.setInputValue("useAnalyticalBDSolution", "false");
        double adb = model.calculateTreeLogLikelihood(tree);
        System.out.println("ADB logL = " + adb); // -643.8421

        assertEquals(adb, bd, tree.getLeafNodeCount() * 0.05);
    }


    // compare log-likelihood curves for the tree and assert close match between implementations
    // see https://github.com/pilarskj/ADB-analysis/accuracy_evaluation
    @Test
    public void testBDMMSeq() throws Exception {

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/treeBD.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // fix all but one parameter
        double C = 5; // lifetime
        double d = 0.1; // death probability
        // double rho = 0.1; // sampling probability
        double origin = 50;

        // loop over different parameter values
        double start = 0.01;
        double end = 1;
        double step = 0.02;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/loglikTreeBD.csv", true);
        DecimalFormat df = new DecimalFormat("0.00");
        for (double i = start; i <= end; i += step) {

            // calculate birth and death rate
            double rho = i;
            double lambda = (1 - d) / C;
            double mu = d / C;

            Parameterization parameterization = new CanonicalParameterization();
            parameterization.initByName(
                    "typeSet", new TypeSet(1),
                    "processLength", new RealParameter(Double.toString(origin)),
                    "birthRate", new SkylineVectorParameter(
                            null,
                            new RealParameter(Double.toString(lambda))),
                    "deathRate", new SkylineVectorParameter(
                            null,
                            new RealParameter(Double.toString(mu))),
                    "birthRateAmongDemes", new SkylineMatrixParameter(null, null),
                    "migrationRate", new SkylineMatrixParameter(null, null),
                    "samplingRate", new SkylineVectorParameter(
                            null,
                            new RealParameter("0")),
                    "rhoSampling", new TimedParameter(
                            new RealParameter(Double.toString(origin)),
                            new RealParameter(Double.toString(rho))),
                    "removalProb", new SkylineVectorParameter(
                            null,
                            new RealParameter("0")));

            BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
            density.initByName(
                    "parameterization", parameterization,
                    "startTypePriorProbs", new RealParameter("1"),
                    "tree", tree,
                    "typeLabel", "type",
                    "parallelize", false,
                    "useAnalyticalSingleTypeSolution", true);

            double logL = density.calculateLogP();
            writer.write("rho," + df.format(i) + "," + logL + ",bdmm\n"); // adapt strings
        }
        writer.close();
    }


    @Test
    public void testADBSeq() throws Exception {

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/treeBD.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // either approx or exact calculations
        boolean approx = true;

        // fix all but one parameter
        double C = 5; // lifetime
        double d = 0.1; // death probability
        // double rho = 0.1; // sampling probability
        double origin = 50;

        // loop over different parameter values
        double start = 0.01;
        double end = 1;
        double step = 0.02;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/ADB-analysis/accuracy_evaluation/loglikTreeBD.csv", true);
        DecimalFormat df = new DecimalFormat("0.00");
        for (double i = start; i <= end; i += step) {
            double rho = i;

            GammaBranchingModel model = new GammaBranchingModel();
            model.initByName(
                    "tree", tree,
                    "shapeInteger", new IntegerParameter("1"),
                    "lifetime", new RealParameter(Double.toString(C)),
                    "deathprob", new RealParameter(Double.toString(d)),
                    "rho", new RealParameter(Double.toString(rho)),
                    "origin", new RealParameter(Double.toString(origin)),
                    "approx", approx,
                    "useAnalyticalBDSolution", false);

            double logL = model.calculateTreeLogLikelihood(tree);
            writer.write("rho," + df.format(i) + "," + logL + ",adb_approx\n"); // adapt strings
        }
        writer.close();
    }
}

