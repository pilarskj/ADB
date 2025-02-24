
package test.adbp;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import bdsky.evolution.speciation.BirthDeathSkylineModel;
import beast.base.evolution.speciation.BirthDeathGernhard08Model;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.text.DecimalFormat;

public class BDTest {

    @Test
    public void testBD() throws Exception {
        // initialize
        BirthDeathGernhard08Model model = new BirthDeathGernhard08Model();

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("type", "labeled");
        model.setInputValue("birthDiffRate", new RealParameter("0.5"));
        model.setInputValue("relativeDeathRate", new RealParameter("0"));
        model.setInputValue("sampleProbability", new RealParameter("0.9"));
        model.setInputValue("originHeight", new RealParameter("15"));

        model.initAndValidate();

        // calculate tree LL
        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL); // value does not match with ADBP (maybe another conditioning?)
    }

    @Test
    public void testBDSKY() throws Exception {
        // initialize
        BirthDeathSkylineModel model = new BirthDeathSkylineModel();

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
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
        System.out.println(logL); // value matches with ADBP without tree factor (has been discarded)
    }


    @Test
    public void testBDMM() throws Exception {

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/calculations/tree_bd.newick",
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
        System.out.println(logL); // value matches with ADBP
    }

    @Test
    public void testBDMMSeq() throws Exception {

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/test2/rhoLow/trees_bd/tree_1.newick",//"/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/calculations/tree_bd.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        /*
        double d = 0.1; // deathprob
        // loop over different parameter values (lifetime)
        double start = 1;
        double end = 25;
        double step = 0.5;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/calculations/loglik_bd.csv", true);
        */
        // loop over different parameter values
        double d = 0.1;
        double start = 4;
        double end = 8;
        double step = 0.05;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/validation/test2/rhoLow/trees_bd/loglik_cf_lifetime.csv", true);
        DecimalFormat df = new DecimalFormat("0.00");
        for (double i = start; i <= end; i += step) {
            double C = i;
            double lambda = (1 - d) / C;
            double mu = d / C;

            // initialize
            Parameterization parameterization = new CanonicalParameterization();
            parameterization.initByName(
                    "typeSet", new TypeSet(1),
                    "processLength", new RealParameter("75.7247271495057"), //50
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
                            new RealParameter("75.7247271495057"), //50
                            new RealParameter("0.001")), //0.1
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
            writer.write(df.format(i) + "," + logL + ",bdmm\n");
        }
        writer.close();


    }
}

