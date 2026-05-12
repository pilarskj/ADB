package adb.operators;

import adb.tree.AnnotatedTreeParser;
import org.junit.jupiter.api.Test;

public class AnnotatedTreeOperatorTest {

    @Test
    public void testWilsonBalding() {

        // initialize tree
        String newick = "(((((((26:0.704651623):0.9815616586):0.9724532301):1.076834454):1.249319317,((((15:0.7965609946):1.052640805):0.9526174057):1.110375181,(((17:0.7945840806):1.02083498):1.106252391):0.9905229342):1.072625896):0.904531321,(((((18:0.7883689605):1.024284023):0.9996663794,((54:0.9184522562):0.9222398679):0.9716272388):0.9727550889):1.034266685):1.070010466):1.036790172):0.9334925396;";
        AnnotatedTreeParser tree = new AnnotatedTreeParser();
        tree.initByName("newick", newick);

        // initialize operator
        AnnotatedWilsonBalding operator = new AnnotatedWilsonBalding();
        operator.initByName("tree", tree, "weight", 1.0);

        // trigger proposal
        operator.proposal();
    }

}
