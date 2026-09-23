package adb;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
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

import java.io.FileWriter;
import java.text.DecimalFormat;


// Compare log-likelihood curves for the tree and assert close match between implementations
// see https://github.com/pilarskj/ADB-analysis/accuracy_evaluation
public class BDCurves {

    void main() throws Exception {

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "src/test/resources/treeBD.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double origin = 50;

        // fix all but one parameter
        double l = 5; // lifetime
        double d = 0.1; // death probability
        // double rho = 0.1; // sampling probability

        // loop over different parameter values
        double start = 0.01;
        double end = 1;
        double step = 0.02;
        FileWriter writer = new FileWriter("src/test/out/loglikTreeBD.csv", true);
        DecimalFormat df = new DecimalFormat("0.00");

        double logL;
        for (double i = start; i <= end; i += step) {

            // calculate birth and death rate
            double rho = i;
            double lambda = (1 - d) / l;
            double mu = d / l;

            // BDMM-Prime likelihood
            Parameterization parameterization = new CanonicalParameterization();
            parameterization.initByName(
                    "typeSet", new TypeSet(1),
                    "processLength", new RealScalarParam<>(origin, NonNegativeReal.INSTANCE),
                    "birthRate", new SkylineVectorParameter(null,
                            new RealVectorParam<>(new double[] {lambda}, NonNegativeReal.INSTANCE), 1),
                    "deathRate", new SkylineVectorParameter(null,
                            new RealVectorParam<>(new double[] {mu}, NonNegativeReal.INSTANCE), 1),
                    "birthRateAmongDemes", new SkylineMatrixParameter(null, null),
                    "migrationRate", new SkylineMatrixParameter(null, null),
                    "samplingRate", new SkylineVectorParameter(null,
                            new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 1),
                    "rhoSampling", new TimedParameter(
                            new RealVectorParam<>(new double[] {origin}, NonNegativeReal.INSTANCE),
                            new RealVectorParam<>(new double[] {rho}, UnitInterval.INSTANCE)),
                    "removalProb", new SkylineVectorParameter(null,
                            new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 1)
            );

            BirthDeathMigrationDistribution bdmm = new BirthDeathMigrationDistribution();
            bdmm.initByName(
                    "tree", tree,
                    "parameterization", parameterization,
                    "startTypePriorProbs",  new SimplexParam(new double[] {1.0}),
                    "useAnalyticalSingleTypeSolution", true);

            logL = bdmm.calculateLogP();
            writer.write("rho," + df.format(i) + "," + logL + ",bdmm\n"); // adapt strings

            // ADB likelihood
            RealScalar<PositiveReal> lifetimeParam = new RealScalarParam<>(l, PositiveReal.INSTANCE);
            RealScalar<PositiveReal> shapeRealParam = new RealScalarParam<>(1.0, PositiveReal.INSTANCE);
            IntScalar<PositiveInt> shapeIntParam = new IntScalarParam<>(1, PositiveInt.INSTANCE);
            RealScalar<UnitInterval> deathprobParam = new RealScalarParam<>(d, UnitInterval.INSTANCE);
            RealScalar<UnitInterval> rhoParam = new RealScalarParam<>(rho, UnitInterval.INSTANCE);

            GammaBranchingModel adbE = new GammaBranchingModel();
            adbE.initByName("tree", tree, "origin", origin,
                    "lifetime", lifetimeParam, "shapeReal", shapeRealParam, "deathprob", deathprobParam, "rho", rhoParam,
                    "approx", false, "useAnalyticalBDSolution", false);

            logL = adbE.calculateTreeLogLikelihood(tree);
            writer.write("rho," + df.format(i) + "," + logL + ",adb_exact\n"); // adapt strings

            GammaBranchingModel adbA = new GammaBranchingModel();
            adbA.initByName("tree", tree, "origin", origin,
                    "lifetime", lifetimeParam, "shapeInteger", shapeIntParam, "deathprob", deathprobParam, "rho", rhoParam,
                    "approx", true, "useAnalyticalBDSolution", false);

            logL = adbA.calculateTreeLogLikelihood(tree);
            writer.write("rho," + df.format(i) + "," + logL + ",adb_approx\n"); // adapt strings
        }

        writer.close();
    }

}
