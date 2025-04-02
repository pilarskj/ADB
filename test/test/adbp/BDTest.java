
package test.adbp;

import adbp.GammaBranchingModel;
import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import bdsky.evolution.speciation.BirthDeathSkylineModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.text.DecimalFormat;

import static org.junit.jupiter.api.Assertions.assertEquals;

// compare estimates between ADBP, BDSKY and BDMM-Prime for shape = 1 (BD case)
public class BDTest {

    @Test
    public void testBDSKY() throws Exception {
        // initialize
        BirthDeathSkylineModel model = new BirthDeathSkylineModel();

        // get tree with 100 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree_bd.newick",
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

        model.initAndValidate();

        // calculate tree LL
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL); // -348.3790 value matches with ADBP without tree factor (has been discarded)
    }


    @Test
    public void testBDMM() throws Exception {

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree_bd.newick",
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
        tree.initByName("fileName", "examples/tree_bd.newick",
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

        assertEquals(adb, bd, 0.5);
    }


    // TO REMOVE (function for comparing likelihood curves)
    @Test
    public void testBDMMSeq() throws Exception {

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/calculations/tree_bd_low_sampling.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double C = 5; // lifetime
        //double d = 0.1; // death probability
        double rho = 0.001; // sampling probability
        // loop over different parameter values
        double start = 0;
        double end = 0.5;
        double step = 0.02;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/calculations/loglik_bd_low_sampling.csv", true);
        DecimalFormat df = new DecimalFormat("0.00");
        for (double i = start; i <= end; i += step) {
            double d = i;
            double lambda = (1 - d) / C;
            double mu = d / C;

            // initialize
            Parameterization parameterization = new CanonicalParameterization();
            parameterization.initByName(
                    "typeSet", new TypeSet(1),
                    "processLength", new RealParameter("100"),
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
                            new RealParameter("100"),
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
            writer.write("d," + df.format(i) + "," + logL + ",bdmm\n");
        }
        writer.close();
    }
}

