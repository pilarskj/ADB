package adb.tree;

import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.Test;

public class AnnotatedTreeTest {

    @Test
    public void testTreeClass() {

        // define tree
        // not typed:
        // String newick = "((D:5.0,C:5.0):6.0,(A:8.0,B:8.0):3.0):0.0;";
        // double origin = 15;

        // tip-typed:
        String newick = "((1[&type=0]:1.378074977,11[&type=1]:1.378074977):2.6813596,(3[&type=0]:1.71378645,6[&type=0]:1.71378645):2.345648127):3.940565423;";
        // origin = 8

        // fully-typed:
        //String newick = "((((1[&type=0]:1.378074977,(11[&type=1]:0.9265429956)[&type=1]:0.4515319812)[&type=1]:0.5929298741)[&type=0]:2.088429726,(((3[&type=0]:0.4513672607)[&type=1]:1.262419189,((6[&type=0]:0.05632331623)[&type=1]:0.7255498933)[&type=1]:0.9319132401)[&type=1]:0.2056271944)[&type=0]:2.140020933)[&type=0]:2.073580344)[&type=0]:1.866985079;";
        // origin = 8

        // parse tree
        AnnotatedTreeParser parser = new AnnotatedTreeParser();
        // TreeParser parser = new TreeParser();
        parser.initByName("newick", newick);
        //parser.initByName("newick",newick, "adjustTipHeights", true, "IsLabelledNewick", true);

        System.out.print(parser);
    }
}