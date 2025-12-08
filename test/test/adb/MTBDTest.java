package test.adb;

import adb.BranchList;
import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import feast.fileio.TreeFromNewickFile;
import org.junit.jupiter.api.Test;

import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.stream.Collectors;

import static adb.MTLogLikelihood.calcMTLogLikelihood;

public class MTBDTest {

    @Test
    public void testMTADB() throws Exception {

        // simulated tree (typed)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/mtADB/trees/tiptree_bdmm.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double t_or = 20;
        int type_or = 0;
        double[] a = {2, 5};
        double[] b = {1, 1};
        double[] d = {0.1, 0.2};
        double rho = 1;
        double[][] Xsi_s = {{0.2, 0}, {0, 0.6}};
        double[][] Xsi_as = {{0, 0.8}, {0.4, 0}};

        // options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 12);
        int mB = (int)Math.pow(2, 10);

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
    public void testMTADBBug() throws Exception {

        // simulated tree (typed)
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "test_data/treeBD.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        double t_or = 50;
        int type_or = 0;
        double[] a = {5};
        double[] b = {1};
        double[] d = {0.1};
        double rho = 0.1;
        double[][] Xsi_s = {{1}};
        double[][] Xsi_as = {{0}};

        // options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);

        // get branches
        BranchList branches = new BranchList(tree, t_or);

        double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_s, Xsi_as, t_or, type_or, branches,
                maxIt, tolP, tolB, mP, mB);

        System.out.println(logL);

        /*// loop over different parameter values
        double start = 1;
        double end = 20;
        double step = 0.5;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/P1_AgeDependentTrees/multi-type/loglik_tree_bd.csv", true);
        DecimalFormat df = new DecimalFormat("0.0");
        for (double x = start; x <= end; x += step) {

            // change theta1
            a[0] = x;

            double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_s, Xsi_as, t_or, type_or, branches,
                    maxIt, tolP, tolB, mP, mB);
            writer.write("lifetime," + df.format(x) + "," + logL + ",mtadb\n"); // adapt strings
        }
        writer.close(); */
    }


    @Test
    public void testBDMM() throws Exception {

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/mtADB/profiling/tree_50.newick", //"/Users/jpilarski/Projects/mtADB/trees/tiptree_bdmm.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // initialize
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", new RealParameter("40"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.09 0.096"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.05 0.04"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.36 0.064"), 2),
                "migrationRate", new SkylineMatrixParameter(null, new RealParameter("0"), 2),
                "samplingRate", new SkylineVectorParameter(null, new RealParameter("0"), 2),
                "rhoSampling", new TimedParameter(
                        new RealParameter("40"),
                        new RealParameter("0.5"), 2),
                "removalProb", new SkylineVectorParameter(null, new RealParameter("1"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", new RealParameter("1 0"),
                "tree", tree,
                "typeLabel", "type",
                "parallelize", true);

        double logL = density.calculateLogP();
        System.out.println(logL);
    }


    @Test
    public void testBDMMSeq() throws Exception {

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/mtADB/trees/tiptree_bdmm.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // fix all but one parameter
        int ntypes = 2;
        double t_or = 20; // 10
        int type_or = 0;
        double[] a = {2, 5};
        double[] d = {0.1, 0.2};
        double rho = 1;
        double[][] Xsi_s = {{0.2, 0}, {0, 0.6}};
        double[][] Xsi_as = {{0, 0.8}, {0.4, 0}};

        // loop over different parameter values
        double start = 0.5;
        double end = 10;
        double step = 0.5;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/mtADB/loglik/loglik_tree_bdmm.csv", false);
        DecimalFormat df = new DecimalFormat("0.0");
        writer.write("parameter,value,logL,method\n");
        for (double x = start; x <= end; x += step) {

            // change parameter
            a[0] = x;

            // calculate BDMM parameters
            double[] h = new double[ntypes];
            h[type_or] = 1;
            double[][] Xsi = new double[ntypes][ntypes];
            for (int i = 0; i < ntypes; i++) {
                for (int j = 0; j < ntypes; j++) {
                    if (i == j) { Xsi[i][j] = Xsi_s[i][j]; }
                    else { Xsi[i][j] = Xsi_as[i][j]; }
                }
            }
            double[] birthRate = new double[ntypes];
            double[][] birthRateAmongDemes = new double[ntypes][ntypes];
            double[] deathRate = new double[ntypes];
            for (int i = 0; i < ntypes; i++) {
                birthRate[i] = Xsi[i][i] * (1 - d[i]) / a[i];
                deathRate[i] = d[i] / a[i];
                for (int j = 0; j < ntypes; j++) {
                    if (i != j) { birthRateAmongDemes[i][j] = Xsi[i][j] * (1 - d[i]) / a[i]; }
                }
            }

            String hS = Arrays.stream(h).mapToObj(Double::toString).collect(Collectors.joining(" "));
            String birthRateS = Arrays.stream(birthRate).mapToObj(Double::toString).collect(Collectors.joining(" "));
            String birthRateAmongDemesS = Arrays.stream(birthRateAmongDemes)
                    .flatMapToDouble(Arrays::stream)
                    .filter(v -> v != 0.0)
                    .mapToObj(Double::toString)
                    .collect(Collectors.joining(" "));
            String deathRateS = Arrays.stream(deathRate).mapToObj(Double::toString).collect(Collectors.joining(" "));

            Parameterization parameterization = new CanonicalParameterization();
            parameterization.initByName(
                    "typeSet", new TypeSet(ntypes),
                    "processLength", new RealParameter(Double.toString(t_or)),
                    "birthRate", new SkylineVectorParameter(
                            null,
                            new RealParameter(birthRateS), ntypes),
                    "deathRate", new SkylineVectorParameter(
                            null,
                            new RealParameter(deathRateS), ntypes),
                    "birthRateAmongDemes", new SkylineMatrixParameter(
                            null,
                            new RealParameter(birthRateAmongDemesS), ntypes),
                    "migrationRate", new SkylineMatrixParameter(null, new RealParameter("0"), ntypes),
                    "samplingRate", new SkylineVectorParameter(null, new RealParameter("0"), ntypes),
                    "rhoSampling", new TimedParameter(
                            new RealParameter(Double.toString(t_or)),
                            new RealParameter(Double.toString(rho)), ntypes),
                    "removalProb", new SkylineVectorParameter(null, new RealParameter("1"), ntypes));

            BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
            density.initByName(
                    "parameterization", parameterization,
                    "startTypePriorProbs", new RealParameter(hS),
                    "tree", tree,
                    "typeLabel", "type",
                    "parallelize", false);

            double logL = density.calculateLogP();
            writer.write("lifetime0," + df.format(x) + "," + logL + ",bdmm\n"); // adapt strings
        }
        writer.close();
    }


    @Test
    public void testMTADBSeq() throws Exception {

        // get tree
        Tree tree = new TreeFromNewickFile();
        tree.initByName("fileName", "/Users/jpilarski/Projects/mtADB/trees/tiptree_bdmm.newick",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        // fix all but one parameter
        int ntypes = 2;
        double t_or = 20;
        int type_or = 0;
        double[] a = {2, 5};
        double[] b = {1, 1};
        double[] d = {0.1, 0.2};
        double rho = 1;
        double[][] Xsi_s = {{0.2, 0}, {0, 0.6}};
        double[][] Xsi_as = {{0, 0.8}, {0.4, 0}};

        // options
        int maxIt = 100;
        double tolP = 1e-12;
        double tolB = 1e-6;
        int mP = (int)Math.pow(2, 14);
        int mB = (int)Math.pow(2, 12);

        // get branches
        BranchList branches = new BranchList(tree, t_or);

        // loop over different parameter values
        double start = 0.5;
        double end = 10;
        double step = 0.5;
        FileWriter writer = new FileWriter("/Users/jpilarski/Projects/mtADB/loglik/loglik_tree_bdmm.csv", true);
        DecimalFormat df = new DecimalFormat("0.0");
        for (double x = start; x <= end; x += step) {

            // change parameter
            a[0] = x;

            double logL = calcMTLogLikelihood(a, b, d, rho, Xsi_s, Xsi_as, t_or, type_or, branches,
                    maxIt, tolP, tolB, mP, mB);
            writer.write("lifetime0," + df.format(x) + "," + logL + ",mtadb\n"); // adapt strings
        }
        writer.close();
    }
}




