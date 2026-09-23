package adb;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.*;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;


// Test likelihood calculation under ADB
public class GammaBranchingModelTest {

    @Test
    public void testGammaBranchingModel() throws Exception {

        // define tree (very small)
        Tree tree = new TreeParser("((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;", false);

        // define parameters
        RealScalar<PositiveReal> lifetime = new RealScalarParam<>(5.0, PositiveReal.INSTANCE);
        RealScalar<UnitInterval> deathprob = new RealScalarParam<>(0.1, UnitInterval.INSTANCE);
        RealScalar<UnitInterval> rho = new RealScalarParam<>(0.1, UnitInterval.INSTANCE);
        double origin = 15.0;

        // define models
        GammaBranchingModel exact = new GammaBranchingModel();
        exact.initByName("tree", tree,
                "lifetime", lifetime, "deathprob", deathprob, "rho", rho,
                "shapeReal", new RealScalarParam<>(5.0, PositiveReal.INSTANCE),
                "origin", origin, "approx", false
        );

        GammaBranchingModel approx = new GammaBranchingModel();
        approx.initByName("tree", tree,
                "lifetime", lifetime, "deathprob", deathprob, "rho", rho,
                "shapeInteger", new IntScalarParam<>(5, PositiveInt.INSTANCE),
                "origin", origin, "approx", true
        );

        // calculate tree log-likelihood
        double likE = exact.calculateTreeLogLikelihood(tree);
        System.out.println("ADB exact logL = " + likE);

        double likA = approx.calculateTreeLogLikelihood(tree);
        System.out.println("ADB approx logL = " + likA);

        assertEquals(likE, likA,  0.05);
    }
}