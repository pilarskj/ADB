package test.adb;

import adb.BranchList;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;

import static adb.MTLogLikColored.calcMTLogLikColored;

public class MTProfilingColored {

    // profile runtime wrt. tree size
    @Test
    public void profileRuntimeSize() throws Exception {

        // define parameters
        // ADB
        int originType = 0;
        double[] lifetime = {2, 5};
        double[] shape = {1, 1};
        double[] death = {0.1, 0.2};
        double rho = 0.5;
        double[][] sTransitions = {{0.2, 0}, {0, 0.6}};
        double[][] asTransitions = {{0, 0.8}, {0.4, 0}};

        int maxIt = 100;
        double tol = 1e-12;
        int m12 = (int)Math.pow(2, 12);
        int m14 = (int)Math.pow(2, 14);
        int m16 = (int)Math.pow(2, 16);

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/mtADB/profiling/colored/tree_data_size.tsv"));
        BufferedWriter bwTime = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/colored/runtimes_size.csv"));
        BufferedWriter bwLik = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/colored/loglik_size.csv"));
        br.readLine(); // skip header in the input
        bwTime.write("tree,stepSize12,stepSize14,stepSize16\n");  // add header in the output
        bwLik.write("tree,stepSize12,stepSize14,stepSize16\n");  // add header in the output

        // read each line from the csv
        String line;
        while ((line = br.readLine()) != null) {
            String[] values = line.split("\t"); // split the line by seperator

            // read values from columns
            int treeNr = Integer.parseInt(values[0]);
            double originTime = Double.parseDouble(values[1]);
            System.out.println(treeNr + ", " + originTime);

            String newick = values[3];
            Tree tree = new TreeParser(newick, true);
            BranchList branches = new BranchList(tree, originTime);

            // record calculation times
            long startTime;
            long endTime;

            startTime = System.nanoTime();
            double ss12Lik = calcMTLogLikColored(lifetime, shape, death, rho, sTransitions, asTransitions, originTime, originType, branches,
                    maxIt, tol, m12);
            endTime = System.nanoTime();
            double ss12Time = (endTime - startTime) / 1e+6;  // divide to ms

            startTime = System.nanoTime();
            double ss14Lik = calcMTLogLikColored(lifetime, shape, death, rho, sTransitions, asTransitions, originTime, originType, branches,
                    maxIt, tol, m14);
            endTime = System.nanoTime();
            double ss14Time = (endTime - startTime) / 1e+6;  // divide to ms

            startTime = System.nanoTime();
            double ss16Lik = calcMTLogLikColored(lifetime, shape, death, rho, sTransitions, asTransitions, originTime, originType, branches,
                    maxIt, tol, m16);
            endTime = System.nanoTime();
            double ss16Time = (endTime - startTime) / 1e+6;  // divide to ms

            bwTime.write(treeNr + "," + ss12Time + "," + ss14Time + "," + ss16Time);
            bwTime.newLine();
            bwLik.write(treeNr + "," + ss12Lik + "," + ss14Lik + "," + ss16Lik);
            bwLik.newLine();
        }
        bwTime.close();
        bwLik.close();
    }


    // profile runtime wrt. number of types
    @Test
    public void profileRuntimeTypes() throws Exception {

        int[] ntypes = new int[15];
        for (int i = 0; i < ntypes.length; i++) {
            ntypes[i] = 2 + i;
        }
        System.out.println(Arrays.toString(ntypes));

        // define parameters
        double rho = 0.5;
        int originType = 0;
        double r = 0.2; // probability of symmetric transition
        int m12 = (int)Math.pow(2, 12);
        int m14 = (int)Math.pow(2, 14);
        int m16 = (int)Math.pow(2, 16);
        int maxIt = 100;
        double tol = 1e-12;

        HashMap<Integer, double[]> lifetime = new HashMap<>();
        HashMap<Integer, double[]> shape = new HashMap<>();
        HashMap<Integer, double[]> death = new HashMap<>();
        HashMap<Integer, double[][]> sTransitions = new HashMap<>();
        HashMap<Integer, double[][]> asTransitions = new HashMap<>();

        for (int n : ntypes) {
            double[] vec = new double[n];
            Arrays.fill(vec, 1.0);
            lifetime.put(n, vec);
            shape.put(n, vec);
            Arrays.fill(vec, 0.1);
            death.put(n, vec);
            double[][] smat = new double[n][n];
            double[][] asmat = new double[n][n];
            for (int i = 0; i < n; i++) {
                smat[i][i] = r; // set diagonal entries
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        asmat[i][j] = (1-r)/(n-1);
                    }
                }
            }
            sTransitions.put(n, smat);
            asTransitions.put(n, asmat);
        }

        // declare in- and out-files
        BufferedReader br = new BufferedReader(new FileReader("/Users/jpilarski/Projects/mtADB/profiling/colored/tree_data_ntypes.tsv"));
        BufferedWriter bwTime = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/colored/runtimes_ntypes.csv"));
        BufferedWriter bwLik = new BufferedWriter(new FileWriter("/Users/jpilarski/Projects/mtADB/profiling/colored/loglik_ntypes.csv"));
        br.readLine(); // skip header in the input
        bwTime.write("tree,stepSize12,stepSize14,stepSize16\n");  // add header in the output
        bwLik.write("tree,stepSize12,stepSize14,stepSize16\n");  // add header in the output

        // read each line from the csv
        String line;
        //int lineNr = 0;
        //int stopLine = 5;
        while ((line = br.readLine()) != null) { // && lineNr < stopLine) {
            //lineNr++;
            String[] values = line.split("\t"); // split the line by seperator

            // read values from columns
            int treeNr = Integer.parseInt(values[0]);
            int n = Integer.parseInt(values[1]);
            double originTime = Double.parseDouble(values[2]);
            System.out.println(treeNr + ", " + n + ", " + originTime);

            String newick = values[3];
            Tree tree = new TreeParser(newick, true);
            BranchList branches = new BranchList(tree, originTime);

            // record calculation times
            long startTime;
            long endTime;

            startTime = System.nanoTime();
            double ss12Lik = calcMTLogLikColored(lifetime.get(n), shape.get(n), death.get(n), rho,
                    sTransitions.get(n), asTransitions.get(n), originTime, originType, branches,
                    maxIt, tol, m12);
            endTime = System.nanoTime();
            double ss12Time = (endTime - startTime) / 1e+6;  // divide to ms

            startTime = System.nanoTime();
            double ss14Lik = calcMTLogLikColored(lifetime.get(n), shape.get(n), death.get(n), rho,
                    sTransitions.get(n), asTransitions.get(n), originTime, originType, branches,
                    maxIt, tol, m14);
            endTime = System.nanoTime();
            double ss14Time = (endTime - startTime) / 1e+6;  // divide to ms

            startTime = System.nanoTime();
            double ss16Lik = calcMTLogLikColored(lifetime.get(n), shape.get(n), death.get(n), rho,
                    sTransitions.get(n), asTransitions.get(n), originTime, originType, branches,
                    maxIt, tol, m16);
            endTime = System.nanoTime();
            double ss16Time = (endTime - startTime) / 1e+6;  // divide to ms

            bwTime.write(treeNr + "," + ss12Time + "," + ss14Time + "," + ss16Time);
            bwTime.newLine();
            bwLik.write(treeNr + "," + ss12Lik + "," + ss14Lik + "," + ss16Lik);
            bwLik.newLine();
        }
        bwTime.close();
        bwLik.close();
    }
}
