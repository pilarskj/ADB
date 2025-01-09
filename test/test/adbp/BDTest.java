package test.adbp;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.speciation.BirthDeathGernhard08Model;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

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
    public void testBDMM() throws Exception {

        // get tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // initialize
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new RealParameter("15"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.4")), //0.5
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1")), //0
                "birthRateAmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("15"),
                        new RealParameter("0.9")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0")));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1"), //startTypePriorProbs - how to use module instead of local installation?
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", true);

        double logL = density.calculateLogP();
        System.out.println(logL); // value matches with ADBP
    }
}
