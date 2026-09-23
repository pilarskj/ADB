package adb;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import bdmmprime.util.ProcessLength;
import beast.base.evolution.tree.Tree;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveInt;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.type.IntScalar;
import beast.base.spec.type.RealScalar;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


// Comparison of likelihood calculation in ADB and BDMM-Prime for tree with shape = 1 (BD case)
public class BDTest {

    @Test
    public void testADBvsBDMM() throws Exception {

        // get tree with 100 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "src/test/resources/treeBD.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // ADB parameters
        double origin = 50.0;
        double lifetime = 5.0;
        double deathprob = 0.1;
        double rho = 0.1;

        // BD parameters
        double birthRate = (1 - deathprob) / lifetime;
        double deathRate = deathprob / lifetime;


        // reference
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new ProcessLength(tree),
                // new RealScalarParam<>(origin, NonNegativeReal.INSTANCE), if not conditionOnRoot
                "birthRate", new SkylineVectorParameter(null,
                        new RealVectorParam<>(new double[] {birthRate}, NonNegativeReal.INSTANCE), 1),
                "deathRate", new SkylineVectorParameter(null,
                        new RealVectorParam<>(new double[] {deathRate}, NonNegativeReal.INSTANCE), 1),
                "birthRateAmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingRate", new SkylineVectorParameter(null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 1),
                "rhoSampling", new TimedParameter(
                        // new RealVectorParam<>(new double[] {origin}, NonNegativeReal.INSTANCE), if not conditionOnRoot
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {rho}, UnitInterval.INSTANCE),
                        new ProcessLength(tree)),
                "removalProb", new SkylineVectorParameter(null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 1)
        );

        BirthDeathMigrationDistribution bdmm = new BirthDeathMigrationDistribution();
        bdmm.initByName(
                "tree", tree,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnRoot", true,
                "useAnalyticalSingleTypeSolution", true);

        double ref = bdmm.calculateLogP();
        System.out.println("BDMM-Prime logL = " + ref); // -643.4968 value matches with ADB (slight difference when conditioning on the root - other definition?)


        // ADB calculations
        double value;
        RealScalar<PositiveReal> lifetimeParam = new RealScalarParam<>(lifetime, PositiveReal.INSTANCE);
        RealScalar<PositiveReal> shapeRealParam = new RealScalarParam<>(1.0, PositiveReal.INSTANCE);
        IntScalar<PositiveInt> shapeIntParam = new IntScalarParam<>(1, PositiveInt.INSTANCE);
        RealScalar<UnitInterval> deathprobParam = new RealScalarParam<>(deathprob, UnitInterval.INSTANCE);
        RealScalar<UnitInterval> rhoParam = new RealScalarParam<>(rho, UnitInterval.INSTANCE);

        GammaBranchingModel bd = new GammaBranchingModel();
        bd.initByName("tree", tree, "origin", origin,
                "lifetime", lifetimeParam, "shapeReal", shapeRealParam, "deathprob", deathprobParam, "rho", rhoParam,
                "approx", false, "useAnalyticalBDSolution", true, "conditionOnRoot", true);

        value = bd.calculateTreeLogLikelihood(tree);
        System.out.println("ADB analytical logL = " + value);
        assertEquals(value, ref,  1e-6);

        GammaBranchingModel adbE = new GammaBranchingModel();
        adbE.initByName("tree", tree, "origin", origin,
                "lifetime", lifetimeParam, "shapeReal", shapeRealParam, "deathprob", deathprobParam, "rho", rhoParam,
                "approx", false, "useAnalyticalBDSolution", false, "conditionOnRoot", true);

        value = adbE.calculateTreeLogLikelihood(tree);
        System.out.println("ADB exact logL = " + value);
        assertEquals(value, ref,  0.5);

        GammaBranchingModel adbA = new GammaBranchingModel();
        adbA.initByName("tree", tree, "origin", origin,
                "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                "approx", true, "useAnalyticalBDSolution", false, "conditionOnRoot", true);

        value = adbA.calculateTreeLogLikelihood(tree);
        System.out.println("ADB approx logL = " + value);
        assertEquals(value, ref,  0.5);
    }

}

