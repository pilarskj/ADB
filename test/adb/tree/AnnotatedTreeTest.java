package adb.tree;

import adb.distribution.ADBTreeDistribution;
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
        //TreeParser parser = new TreeParser();
        //parser.initByName("newick",newick, "adjustTipHeights", true, "IsLabelledNewick", true);
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

        tree = aTree.convertAnnotatedTree(false);
        System.out.println(tree);
    }

    @Test
    public void testInitialisationComplex() throws Exception {

        Tree tree = new TreeParser("((1849[&type=2]:6.934016968,(1392[&type=1]:4.353209345,1808[&type=1]:4.353209345):2.580807623):2.599236762,((931[&type=1]:8.518002011,(1104[&type=2]:7.689335149,(2705[&type=0]:3.622823942,890[&type=2]:3.622823942):4.066511207):0.8286668624):0.4522656605,((2768[&type=0]:4.956776144,2014[&type=0]:4.956776144):3.478486295,((573[&type=1]:7.406632908,(439[&type=2]:5.929580764,1697[&type=3]:5.929580764):1.477052145):0.4893282464,((1215[&type=2]:3.953081826,944[&type=3]:3.953081826):3.406588011,(434[&type=2]:5.293885734,(2129[&type=0]:3.600093916,1253[&type=1]:3.600093916):1.693791818):2.065784103):0.5362913174):0.5393012838):0.5350052332):0.5629860586);",
                true);
        Parameterization model = new Parameterization();
        model.initByName("nTypes", 4,
                "lifetime", new RealParameter("0.5 1 2 3"),
                "shape", new IntegerParameter("100 50 20 5"),
                "death", new RealParameter("0.05 0.1 0.1 0.1"),
                "symTransitions", new RealParameter("0.4 0 0 0 0 0.2 0 0 0 0 1 0 0 0 0 1"),
                "asymTransitions", new RealParameter("0 0.6 0 0 0 0 0.5 0.3 0 0 0 0 0 0 0 0"),
                "sampling", new RealParameter("0.005 0.005 0.005 0.005"),
                "originTime", 10.0);

        AnnotatedTree aTree = new AnnotatedTreeInitialiser();
        aTree.initByName("tree", tree, "parameterization", model, "scale", 0.8);
        System.out.println(aTree);

        tree = aTree.convertAnnotatedTree(true);
        System.out.println(tree);

        // calculate tree log-likelihood (check if not -Infinity)
        ADBTreeDistribution distribution = new ADBTreeDistribution();
        distribution.initByName("tree", aTree,
                "parameterization", model);

        double logL = distribution.calculateTreeLogLikelihood(aTree);
        System.out.println(logL);
    }
}