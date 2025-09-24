package test.adb;

import adb.BranchList;
import beast.base.evolution.tree.Tree;
import feast.fileio.TreeFromNewickFile;

import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.util.Arrays;

import static adb.MTLogLikelihood.calcMTLogLikelihood;

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
        BranchList branches = new BranchList(tree, t_or, type_or);

        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();
        System.out.println(n_int + " internal and " + n_ext + " external branches");

        double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_s, Xsi_as, t_or, type_or, branches,
                maxIt, tolP, tolB, mP, mB);

        System.out.println(logL);
    }


    @Test
    public void testMTLogLikelihoodGrid() throws Exception {

        // simulated tree (typed)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/P1_AgeDependentTrees/multi-type/tiptree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double t_or = 10;
        int type_or = 0;
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
        BranchList branches = new BranchList(tree, t_or, type_or);

        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();
        System.out.println(n_int + " internal and " + n_ext + " external branches");

        // loop over different lifetimes
        double[] a;
        double start = 0.5;
        double end = 5;
        double step = 0.5;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/multi-type/logL_lifetime_3D.csv");
        writer.write("lifetime0,lifetime1,logL\n");
        for (double i = start; i <= end; i += step) {
            for (double j = start; j <= end; j += step) {
                a = new double[]{i/b[0], j/b[1]};
                System.out.println(Arrays.toString(a));

                double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_s, Xsi_as, t_or, type_or, branches,
                        maxIt, tolP, tolB, mP, mB);
                System.out.println(logL);

                writer.write(i + "," + j + "," + logL + "\n");
            }
        }
        writer.close();
    }
}