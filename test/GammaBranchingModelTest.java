import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class GammaBranchingModelTest  {

    @Test
    public void testGammaBranching() throws Exception {

        // initialize
        GammaBranchingModel model = new GammaBranchingModel();

        // define tree
        Tree tree = new TreeParser("((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;", false);
        model.setInputValue("tree", tree);

        // set parameters
        model.setInputValue("scale", new RealParameter("1"));
        model.setInputValue("shape", new RealParameter("1"));
        model.setInputValue("rho", new RealParameter("0.1"));
        model.setInputValue("origin", new RealParameter("12"));

        model.initAndValidate();

        double logL = model.calculateTreeLogLikelihood(tree);
        System.out.println(logL);
        assertEquals(logL, -26.83856, 0.00001);
    }
}