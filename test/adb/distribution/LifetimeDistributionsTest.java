package adb.distribution;

import adb.util.Utils;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;

import org.junit.Before;
import org.junit.Test;

public class LifetimeDistributionsTest {

    RealParameter lifetime;
    Parameter shape;
    double[] timeArray;

    LifetimeDistributions distributions;


    @Before
    public void setUp() {
        lifetime = new RealParameter("1.0 1.0 1.0 1.0");
        shape = new RealParameter("10.0 10.0 10.0 10.0");

        double origin = 10.0;
        int nSteps = (int) Math.pow(2, 14);
        timeArray = Utils.linSpace(0, origin, nSteps);

        // Distribution setup
        distributions = new LifetimeDistributions();
        distributions.initByName("lifetime", lifetime, "shape", shape, "timeArray", timeArray);
    }


    @Test
    public void testManualUpdate() {

        // Register all parameters that the MCMC would move
        //State state = new State();
        //state.initByName("stateNode", lifetime, "stateNode", shape);

        // Back-up current state
        //state.store(0);
        //distributions.store();

        // Simulate an MCMC Move (The "Proposal")
        // MUST call startEditing() on any StateNode before changing it
        lifetime.startEditing(null);
        lifetime.setValue(0, 5.0); // Change the first type's lifetime

        // Trigger Calculation
        distributions.requiresRecalculation(); // should return true
        distributions.getLifetimeDistributions();

        /* // Reject move
        //state.restore();
        lifetime.restore();
        distribution.restore();
        distribution.getLifetimeDistributions(); */

        distributions.accept();

        System.out.println("");
    }

}
