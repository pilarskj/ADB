package test.adbp;

import adbp.Branch;
import adbp.BranchList;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import feast.fileio.TreeFromNewickFile;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static adbp.MTLogLikelihood.calcLogLikelihood;

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
        */

        // tree with 150 tips
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "examples/tree.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // params
        int ntypes = 2;
        int ntips = tree.getLeafNodeCount();

        double[] a = {1, 2};
        double[] b = {1, 2};
        double[] d = {0.1, 0.2};
        double rho = 0.9;
        double[][] Xsi_as = {{0, 0.1},
                             {0.2, 0}};
        double[][] Xsi_s = {{0.5, 0.3},
                            {0.1, 0.5}};
        double t_or = 15;
        //int type_or = 0;

        // options
        int m = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);
        int maxit = 100;

        // get branches
        BranchList branches = new BranchList(tree, t_or);
        branches.assignTypes(ntypes);

        int n_int = branches.countInternalBranches();
        int n_ext = branches.countExternalBranches();

        double[] int_s = new double[n_int];
        double[] int_e = new double[n_int];
        int[] left_child = new int[n_int];
        int[] right_child = new int[n_int];
        double[] ext_e = new double[n_ext];
        int[] types = new int[n_int + n_ext];

        // collect info for internal branches
        for (int i = 0; i < n_int; i++) {
            Branch branch = branches.getBranchByIndex(i);
            int_s[i] = branch.startTime;
            int_e[i] = branch.endTime;
            left_child[i] = branch.leftIndex;
            right_child[i] = branch.rightIndex;
            types[i] = branch.nodeType;
        }

        // collect end times for external branches
        for (int i = 0; i < n_ext; i++) {
            Branch branch = branches.getBranchByIndex(i + n_int);
            ext_e[i] = branch.endTime;
            types[i + n_int] = branch.nodeType;
        }
        /*
        System.out.println(Arrays.toString(int_s));
        System.out.println(Arrays.toString(int_e));
        System.out.println(Arrays.toString(left_child));
        System.out.println(Arrays.toString(right_child));
        System.out.println(Arrays.toString(ext_e));
        */

        double logL = calcLogLikelihood(a, b, d, rho, Xsi_as, Xsi_s, t_or, types, int_s, int_e, ext_e,
                left_child, right_child, m, mB, maxit);

        System.out.println(logL);
    }
}
