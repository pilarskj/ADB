package adb.distribution;

import adb.util.Utils;
import beast.base.inference.parameter.RealParameter;
import org.junit.Before;
import org.junit.Test;

public class P0SystemTest {

    RealParameter lifetime;
    P0System P0;
    LifetimeDistributions lifetimeDistributions;

    @Before
    public void setUp() {

        lifetime = new RealParameter("2 5");
        RealParameter shape = new RealParameter("1 1");
        double origin = 20.0;

        Parameterization parameterization = new Parameterization();
        parameterization.initByName("nTypes", 2,
                "lifetime", lifetime,
                "shape", shape,
                "death", new RealParameter("0.1 0.2"),
                "symTransitions", new RealParameter("0.2 0 0 0.6"),
                "asymTransitions", new RealParameter("0 0.8 0.4 0"),
                "sampling", new RealParameter("0.5 0.5"),
                "originTime", origin,
                "originType", 0);

        int maxIt = 100;
        double tol = 1e-12;
        int nSteps = (int) Math.pow(2, 8);
        double[] timeArray = Utils.linSpace(0, origin, nSteps);
        double timeStep = origin / nSteps;

        lifetimeDistributions = new LifetimeDistributions();
        lifetimeDistributions.initByName("lifetime", lifetime, "shape", shape, "timeArray", timeArray);

        // CalculationNode setup
        P0 = new P0System();
        P0.initByName("parameterization", parameterization,
                "lifetimeDistributions", lifetimeDistributions,
                "maxIterations", maxIt,
                "tolerance", tol,
                "nSteps", nSteps,
                "timeArray", timeArray,
                "timeStep", timeStep);
    }


    @Test
    public void testManualUpdate() {

        lifetimeDistributions.store();
        P0.store();

        // Simulate an MCMC Move (The "Proposal")
        // MUST call startEditing() on any StateNode before changing it
        lifetime.startEditing(null);
        lifetime.setValue(0, 5.0); // Change the first type's lifetime

        // Trigger Calculation
        lifetimeDistributions.requiresRecalculation();
        P0.requiresRecalculation();
        P0.getP0();

        // Reject move
        lifetime.restore();
        lifetimeDistributions.restore();
        P0.restore();

        /* // Accept move
        lifetimeDistributions.accept();
        lifetimeDistributions.store();
        P0.accept();
        P0.store(); */

        System.out.println("");
    }

}
