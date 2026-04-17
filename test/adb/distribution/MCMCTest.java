package adb.distribution;

import adb.util.Utils;
import beast.base.core.Input;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.inference.*;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.List;
import java.util.Random;

/* public class MCMCTest {

    // Tree tree;
    RealParameter lifetime;
    Parameter shape;
    double[] timeArray;

    TransformedDistribution distribution;

    ScaleOperator lifetimeScaler;
    ScaleOperator shapeScaler;


    @Before
    public void setUp() {

        // Create tree
        tree = new TreeParser();
        tree.initByName(
                "newick", "((1849[&type=2]:6.934016968,(1392[&type=1]:4.353209345,1808[&type=1]:4.353209345):2.580807623):2.599236762,((931[&type=1]:8.518002011,(1104[&type=2]:7.689335149,(2705[&type=0]:3.622823942,890[&type=2]:3.622823942):4.066511207):0.8286668624):0.4522656605,((2768[&type=0]:4.956776144,2014[&type=0]:4.956776144):3.478486295,((573[&type=1]:7.406632908,(439[&type=2]:5.929580764,1697[&type=3]:5.929580764):1.477052145):0.4893282464,((1215[&type=2]:3.953081826,944[&type=3]:3.953081826):3.406588011,(434[&type=2]:5.293885734,(2129[&type=0]:3.600093916,1253[&type=1]:3.600093916):1.693791818):2.065784103):0.5362913174):0.5393012838):0.5350052332):0.5629860586);",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // Set parameters
        lifetime = new RealParameter("1 1 1 1");
        shape = new RealParameter("10 10 10 10");

        double origin = 10.0;
        int nSteps = (int) Math.pow(2, 14);
        timeArray = Utils.linSpace(0, origin, nSteps);

        // Create distribution
        System.setProperty("java.only", "true");
        distribution = new TransformedDistribution();
        distribution.initByName(
                "lifetime", lifetime,
                "shape", shape,
                "timeArray", timeArray);

        // Create operators
        lifetimeScaler = new ScaleOperator();
        lifetimeScaler.initByName("parameter", lifetime, "weight", 1.0);
        shapeScaler = new ScaleOperator();
        shapeScaler.initByName("parameter", shape, "weight", 1.0);

        // Initialize the dummy likelihood
        likelihood = new DummyLikelihood();
        likelihood.initByName("dist", distribution);
    }


    @Test
    public void TestMCMC() {
        // Force initialization
        likelihood.calculateLogP();

        // set up state:
        State state = new State();
        state.initByName("stateNode", lifetime, "stateNode", shape);

        // Set up logger:
        Logger lifetimeLogger = new Logger();
        Logger shapeLogger = new Logger();
        lifetimeLogger.initByName("log", lifetime);
        shapeLogger.initByName("log", shape);

        // Set up MCMC:
        MCMC mcmc = new MCMC();
        CompoundDistribution likelihoodForMCMC = new CompoundDistribution();
        likelihoodForMCMC.initByName("distribution", likelihood);
        mcmc.initByName(
                "chainLength", "100",
                "state", state,
                "distribution", likelihoodForMCMC,
                "operator", lifetimeScaler,
                "operator", shapeScaler,
                "logger", lifetimeLogger,
                "logger", shapeLogger);

        // Run MCMC:
        try {
            mcmc.run();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (SAXException e) {
            e.printStackTrace();
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        }

    }
} */
