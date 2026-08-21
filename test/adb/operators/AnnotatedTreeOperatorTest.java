package adb.operators;

import adb.distribution.Parameterization;
import adb.tree.AnnotatedTreeParser;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

public class AnnotatedTreeOperatorTest {

    String newick = "(((((((5:0.704651623):0.9815616586):0.9724532301):1.076834454):1.249319317,((((2:0.7965609946):1.052640805):0.9526174057):1.110375181,(((3:0.7945840806):1.02083498):1.106252391):0.9905229342):1.072625896):0.904531321,(((((4:0.7883689605):1.024284023):0.9996663794,((1:0.9184522562):0.9222398679):0.9716272388):0.9727550889):1.034266685):1.070010466):1.036790172):0.9334925396;";
    AnnotatedTreeParser tree;
    RealParameter lifetime;
    IntegerParameter shape;
    RealParameter death;
    RealParameter sampling;
    Parameterization model;
    Tree flatTree;
    double logHR;

    @BeforeEach
    public void setUp() {
        tree = new AnnotatedTreeParser();
        tree.initByName("newick", newick);

        lifetime = new RealParameter("1.0");
        lifetime.setBounds(0.0, Double.POSITIVE_INFINITY);
        shape = new IntegerParameter("100");
        shape.setBounds(1, 1000);

        death = new RealParameter("0.1");
        death.setBounds(0.0, 1.0);
        sampling = new RealParameter("0.1");
        sampling.setBounds(0.0, 1.0);

        model = new Parameterization();
        model.initByName("nTypes", 1,
                "lifetime", lifetime,
                "shape", shape,
                "death", death, //new RealParameter("0.1"),
                "sampling", sampling, //new RealParameter("0.1"),
                "originTime", 7.859634);
    }

    @AfterEach
    public void printResult() {
        System.out.println(logHR);
        if (logHR != Double.NEGATIVE_INFINITY) {
            flatTree = tree.convertAnnotatedTree(false);
            System.out.println(flatTree.getRoot().toNewick());
        }
    }


    @Test
    public void testEventSampler() {

        // initialize operator
        EventSampler operator = new EventSampler();
        operator.initByName("tree", tree, "parameterization", model, "drawEventCount", true, "proportionBranches", 1.0, "weight", 1.0);

        // trigger proposal
        logHR = operator.proposal();
    }

    @Test
    public void testJointEventSampler() {

        // initialize operator
        JointEventSampler operator = new JointEventSampler();
        //operator.initByName("tree", tree, "parameterization", model, "realParameter", lifetime, "scaleFactor", 0.5, "drawEventCount", true, "weight", 1.0);
        //operator.initByName("tree", tree, "parameterization", model, "intParameter", shape, "windowSize", 50, "drawEventCount", false, "weight", 1.0);
        operator.initByName("tree", tree, "parameterization", model, "parameterInverse", lifetime, "parameterInverse", sampling, "scaleFactor", 0.2, "weight", 1.0);

        // trigger proposal
        logHR = operator.proposal();
    }


    @Test
    public void testWilsonBalding() {

        // initialize operator
        AnnotatedWilsonBalding operator = new AnnotatedWilsonBalding();
        operator.initByName("tree", tree, "parameterization", model, "weight", 1.0);

        // trigger proposal
        logHR = operator.proposal();
    }


    @Test
    public void testExchange() {

        // initialize operator
        AnnotatedExchange operator = new AnnotatedExchange();
        operator.initByName("tree", tree, "parameterization", model, "isNarrow", true, "weight", 1.0);

        // trigger proposal
        logHR = operator.proposal();

        // trigger proposal
        operator.proposal();
    }


    @Test
    public void testSubtreeSlide() {

        // initialize operator
        AnnotatedSubtreeSlide operator = new AnnotatedSubtreeSlide();
        operator.initByName("tree", tree, "parameterization", model, "weight", 1.0, "size", 3);

        // trigger proposal
        logHR = operator.proposal();
    }

    @Test
    public void testScale() {

        // initialize operator
        AnnotatedScale operator = new AnnotatedScale();
        operator.initByName("tree", tree, "parameterization", model, "parameterInverse", lifetime, "weight", 1.0, "scaleFactor", 0.5, "drawEvents", true);

        // trigger proposal
        logHR = operator.proposal();
    }

}
