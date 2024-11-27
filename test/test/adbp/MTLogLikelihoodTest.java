package test.adbp;

import adbp.Branch;
import adbp.BranchList;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;

import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.util.Arrays;

import static adbp.MTLogLikelihood.calcMTLogLikelihood;

public class MTLogLikelihoodTest {

    @Test
    public void testMTLogLikelihood() throws Exception {

        /*
        // small tree
        String newick = "(C:11.0,(A:8.0,B:8.0):3.0):0.0;";
        TreeParser tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        // tree with 150 tips (untyped)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);
         */

        // simulated tree with 26 tips (typed)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree_typed.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // params
        int ntypes = 2;

        double t_or = 15;
        double[] a = {1, 2};
        double[] b = {1, 2};
        double[] d = {0.2, 0.1};
        double rho = 0.5;
        double[][] Xsi_as = {{0, 0.1},
                             {0.2, 0}};
        double[][] Xsi_s = {{0.5, 0.3},
                            {0.1, 0.5}};
        // int type_or = 0;

        // options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);

        // get branches
        BranchList branches = new BranchList(tree, t_or, true);
        // branches.assignRandomTypes(ntypes);

        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();
        // System.out.println(n_int + " internal and " + n_ext + " external branches");

        double[] intS = new double[n_int];
        double[] intE = new double[n_int];
        int[] left_child = new int[n_int];
        int[] right_child = new int[n_int];
        double[] extE = new double[n_ext];
        int[] types = new int[n_int + n_ext];

        // collect info for internal branches
        for (int i = 0; i < n_int; i++) {
            Branch branch = branches.getBranchByIndex(i);
            intS[i] = branch.startTime;
            intE[i] = branch.endTime;
            left_child[i] = branch.leftIndex;
            right_child[i] = branch.rightIndex;
            types[i] = branch.nodeType;
        }

        // collect end times for external branches
        for (int i = 0; i < n_ext; i++) {
            Branch branch = branches.getBranchByIndex(i + n_int);
            extE[i] = branch.endTime;
            types[i + n_int] = branch.nodeType;
        }

        /*
        System.out.println(Arrays.toString(intS));
        System.out.println(Arrays.toString(intE));
        System.out.println(Arrays.toString(left_child));
        System.out.println(Arrays.toString(right_child));
        System.out.println(Arrays.toString(extE));
        */

        double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_as, Xsi_s, t_or, types, intS, intE, extE, left_child, right_child,
                maxIt, tolP, tolB, mP, mB);

        System.out.println(logL);
    }


    @Test
    public void testMTLogLikelihoodSeq() throws Exception {

        // simulated tree with 26 tips (typed)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree_typed.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double t_or = 15;
        double[] a = {1, 2};
        double[] b = {1, 2};
        double[] d = {0.2, 0.1};
        //double rho = 0.5;
        double[][] Xsi_as = {{0, 0.1},
                {0.2, 0}};
        double[][] Xsi_s = {{0.5, 0.3},
                {0.1, 0.5}};

        // options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);

        // get branches
        BranchList branches = new BranchList(tree, t_or, true);
        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();

        double[] intS = new double[n_int];
        double[] intE = new double[n_int];
        int[] left_child = new int[n_int];
        int[] right_child = new int[n_int];
        double[] extE = new double[n_ext];
        int[] types = new int[n_int + n_ext];

        // collect info for internal branches
        for (int i = 0; i < n_int; i++) {
            Branch branch = branches.getBranchByIndex(i);
            intS[i] = branch.startTime;
            intE[i] = branch.endTime;
            left_child[i] = branch.leftIndex;
            right_child[i] = branch.rightIndex;
            types[i] = branch.nodeType;
        }

        // collect end times for external branches
        for (int i = 0; i < n_ext; i++) {
            Branch branch = branches.getBranchByIndex(i + n_int);
            extE[i] = branch.endTime;
            types[i + n_int] = branch.nodeType;
        }

        // loop over different scales
        double start = 0;
        double end = 1;
        double step = 0.1;
        FileWriter writer = new FileWriter("examples/logL_typed_rho.csv");
        writer.write("rho,logL\n");
        for (double i = start; i <= end; i += step) {
            //for (double j = start; j <= end; j += step) {
                double rho = i; //, new double[]{i, j};
                System.out.println(rho); //Arrays.toString(a));
                try {
                    double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_as, Xsi_s, t_or, types, intS, intE, extE, left_child, right_child,
                            maxIt, tolP, tolB, mP, mB);
                    writer.write(i + "," + logL + "\n"); //+ j + ","
                } catch (Exception e) {
                    writer.write(i + "," + "\n"); //+ j + ","
                }
            //}
        }
        writer.close();
    }
}
