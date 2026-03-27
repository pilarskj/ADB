package adb.distribution;

import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.Test;

public class ParameterizationTest {

    // test model parameters and settings
    @Test
    public void testParametrization() throws Exception {
        Parameterization model = new Parameterization();
        model.initByName("nTypes", 2,
                "lifetime", new RealParameter("2 5"),
                "shape", new RealParameter("1 1"),
                "death", new RealParameter("0.1 0.2"),
                "transitions", new RealParameter("0.2 0.8 0.4 0.6"),
                //"sTransitions", new RealParameter("0.2 0 0 0.6"),
                //"asTransitions", new RealParameter("0 0.8 0.4 0"),
                "sampling", new RealParameter("0.5 0.5"),
                "originTime", 20.0,
                "originType", 0);

        System.out.println(model);
    }

}
