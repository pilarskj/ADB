package test.adb;

import adb.BranchList;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import feast.fileio.TreeFromNewickFile;

import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Arrays;

import static adb.MTLogLikelihood.calcMTLogLikelihood;
import static adb.MTLogLikColored.calcMTLogLikColored;

public class MTLogLikelihoodTest {

    @Test
    public void testMTLogLikelihood() throws Exception {

        // simulated tree (typed)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/multi-type/tiptree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double t_or = 10;
        int type_or = 0;
        double[] a = {0.02, 0.2};
        double[] b = {100, 5};
        double[] d = {0.05, 0.1};
        double rho = 0.5;
        double[][] Xsi_s = {{0.5, 0.3}, {0.1, 0.5}};
        double[][] Xsi_as = {{0, 0.1}, {0.2, 0}};

        // options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);

        // get branches
        BranchList branches = new BranchList(tree, t_or);

        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();
        System.out.println(n_int + " internal and " + n_ext + " external branches");

        double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_s, Xsi_as, t_or, type_or, branches,
                maxIt, tolP, tolB, mP, mB);

        System.out.println(logL);
    }


    @Test
    public void testMTLogLikColored() throws Exception {

        // simulated tree (branch-typed)
        String newick = "((((1[&type=0]:1.378074977,(11[&type=1]:0.9265429956)[&type=1]:0.4515319812)[&type=1]:0.5929298741)[&type=0]:2.088429726,(((3[&type=0]:0.4513672607)[&type=1]:1.262419189,((6[&type=0]:0.05632331623)[&type=1]:0.7255498933)[&type=1]:0.9319132401)[&type=1]:0.2056271944)[&type=0]:2.140020933)[&type=0]:2.073580344)[&type=0]:1.866985079;";
        TreeParser tree = new TreeParser();
        tree.initByName("newick", newick, "IsLabelledNewick", true, "adjustTipHeights", true);

        double originTime = 8;
        int originType = 0;

        double[] a = {0.02, 0.2};
        double[] b = {100, 5};
        double[] d = {0.1, 0.2};
        double rho = 0.5;
        double[][] Xsi_s = {{0.3, 0.3}, {0.1, 0.5}};
        double[][] Xsi_as = {{0, 0.4}, {0.4, 0}};

        // options
        int maxIt = 100;
        double tolP = 1e-12;
        int mP = (int)Math.pow(2, 14);

        // get branches
        BranchList branches = new BranchList(tree, originTime);

        double logL = calcMTLogLikColored(a, b, d, rho, Xsi_s, Xsi_as, originTime, originType, branches,
                maxIt, tolP, mP);

        System.out.println(logL);
    }


    @Test
    public void testMTLogLikColoredSize() throws Exception {

        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/mtADB/profiling/colored/tree_49_size.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double originTime = 40;
        int originType = 0;

        double[] a = {2, 5};
        double[] b = {1, 1};
        double[] d = {0.1, 0.2};
        double rho = 0.5;
        double[][] Xsi_s = {{0.2, 0}, {0, 0.6}};
        double[][] Xsi_as = {{0, 0.8}, {0.4, 0}};

        // options
        int maxIt = 100;
        double tol = 1e-12;
        int m = (int)Math.pow(2, 14);

        // get branches
        BranchList branches = new BranchList(tree, originTime);

        double logL = calcMTLogLikColored(a, b, d, rho, Xsi_s, Xsi_as, originTime, originType, branches,
                maxIt, tol, m);

        System.out.println(logL);
    }



    @Test
    public void testMTLogLikelihoodGrid() throws Exception {

        /*// simulated tree (typed)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/multi-type/tiptree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true); */

        // simulated tree (branch-typed)
        String newick = "((((1[&type=0]:1.378074977,(11[&type=1]:0.9265429956)[&type=1]:0.4515319812)[&type=1]:0.5929298741)[&type=0]:2.088429726,(((3[&type=0]:0.4513672607)[&type=1]:1.262419189,((6[&type=0]:0.05632331623)[&type=1]:0.7255498933)[&type=1]:0.9319132401)[&type=1]:0.2056271944)[&type=0]:2.140020933)[&type=0]:2.073580344)[&type=0]:1.866985079;";
        TreeParser tree = new TreeParser();
        tree.initByName("newick", newick, "IsLabelledNewick", true, "adjustTipHeights", true);

        double originTime = 8;
        int originType = 0;

        double[] b = {100, 5};
        double[] d = {0.1, 0.2};
        double rho = 0.5;
        double[][] Xsi_s = {{0.3, 0.3}, {0.1, 0.5}};
        double[][] Xsi_as = {{0, 0.4}, {0.4, 0}};

        // options
        int maxIt = 100;
        double tolP = 1e-12;
        //double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        //int mB = (int)Math.pow(2, 12);

        // get branches
        BranchList branches = new BranchList(tree, originTime);

        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();
        System.out.println(n_int + " internal and " + n_ext + " external branches");

        // loop over different lifetimes
        double[] a;
        double start = 0.8;
        double end = 4;
        double step = 0.2;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/mtADB/loglik/loglik_colored_lifetime_grid.csv");
        DecimalFormat df = new DecimalFormat("0.0");
        writer.write("lifetime0,lifetime1,logL\n");
        for (double i = start; i <= end; i += step) {
            for (double j = start; j <= end; j += step) {
                a = new double[]{i/b[0], j/b[1]};
                System.out.println(Arrays.toString(a));

                // double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_s, Xsi_as, t_or, type_or, branches,
                //         maxIt, tolP, tolB, mP, mB);
                double logL = calcMTLogLikColored(a, b, d, rho, Xsi_s, Xsi_as, originTime, originType, branches,
                        maxIt, tolP, mP);
                System.out.println(logL);

                writer.write(df.format(i) + "," + df.format(j) + "," + logL + "\n");
            }
        }
        writer.close();
    }
}