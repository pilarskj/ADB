package adb.tree;

import adb.distribution.Parameterization;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static adb.util.Utils.distributeDirichlet;

public class AnnotatedTreeTest {

    @Test
    public void testTreeClass() {

        // define tree
        // not typed:
        // String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";

        // tip-typed:
        // String newick = "((1[&type=0]:1.378074977,11[&type=1]:1.378074977):2.6813596,(3[&type=0]:1.71378645,6[&type=0]:1.71378645):2.345648127):3.940565423;";

        // fully-typed:
        String newick = "((((1[&type=0]:1.378074977,(11[&type=1]:0.9265429956)[&type=1]:0.4515319812)[&type=1]:0.5929298741)[&type=0]:2.088429726,(((3[&type=0]:0.4513672607)[&type=1]:1.262419189,((6[&type=0]:0.05632331623)[&type=1]:0.7255498933)[&type=1]:0.9319132401)[&type=1]:0.2056271944)[&type=0]:2.140020933)[&type=0]:2.073580344)[&type=0]:1.866985079;";

        // parse tree
        // TreeParser parser = new TreeParser();
        // parser.initByName("newick",newick, "adjustTipHeights", true, "IsLabelledNewick", true);
        AnnotatedTreeParser parser = new AnnotatedTreeParser();
        parser.initByName("newick", newick);

        System.out.print(parser);
    }

    @Test
    public void testDirichlet() {
        double[] x = {4, 4, 4};
        double alpha = 5;
        double[] res = distributeDirichlet(x, alpha);
        System.out.println(Arrays.toString(res));
    }

    @Test
    public void testInitialisation() throws Exception {

        Tree tree = new TreeParser("((1:43.08956348,17:43.08956348):7.807559093,(8:17.26568023,10:17.26568023):33.63144235):8.415882985;", true);
        Parameterization model = new Parameterization();
        model.initByName("nTypes", 1,
                "lifetime", new RealParameter("10"),
                "shape", new IntegerParameter("20"),
                "death", new RealParameter("0.1"),
                "sampling", new RealParameter("0.1"),
                "originTime", 70.0);

        AnnotatedTree aTree = new AnnotatedTreeInitialiser();
        aTree.initByName("tree", tree, "parameterization", model);
        System.out.println(aTree);
    }
}