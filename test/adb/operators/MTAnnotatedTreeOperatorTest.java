package adb.operators;

import adb.distribution.Parameterization;
import adb.tree.AnnotatedTree;
import adb.tree.AnnotatedTreeParser;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

public class MTAnnotatedTreeOperatorTest {

    String newick = "((((1[&type=3]:0.5050647894)[&type=3]:2.629032343,(((3[&type=3]:0.02658926371)[&type=1]:1.167256789)[&type=1]:0.9582704094)[&type=1]:0.9819806705)[&type=1]:0.8807615117,((((2[&type=1]:0.3540064826)[&type=1]:0.9311262134)[&type=1]:1.069878825)[&type=1]:1.123534175)[&type=0]:0.536312948)[&type=0]:0.518395086,((((((((10[&type=0]:0.2977362855)[&type=0]:0.5473791622)[&type=0]:0.5281869737)[&type=0]:0.4965747618)[&type=0]:0.536755725)[&type=0]:0.4893282464,((((8[&type=2]:0.1301792254)[&type=1]:1.243179893)[&type=0]:0.4946663287)[&type=0]:0.4916443901)[&type=0]:0.5362913174)[&type=0]:0.5393012838)[&type=0]:0.5350052332,(((((5[&type=2]:0.173317248)[&type=2]:1.486935976)[&type=1]:0.8881744641)[&type=0]:0.4923866848,((((7[&type=1]:0.4430713129)[&type=1]:1.124596537)[&type=0]:0.4915929058)[&type=0]:0.5159186473)[&type=0]:0.4656349701)[&type=0]:0.4771876381,((((4[&type=3]:0.6597367559)[&type=1]:0.9506927264)[&type=1]:1.078905666,(((6[&type=1]:0.3101182183)[&type=1]:0.7872153939)[&type=1]:1.103073858)[&type=0]:0.4889276784)[&type=0]:0.4354935528,(((((9[&type=3]:0.1416000553)[&type=1]:0.9612053662,((12[&type=0]:0.0160933408)[&type=0]:0.4995183673)[&type=0]:0.5871937134)[&type=0]:0.5595680455,(((11[&type=1]:0.1244051404)[&type=0]:0.4221686028)[&type=0]:0.5657392679)[&type=0]:0.550060456)[&type=0]:0.4531051483)[&type=0]:0.4907869725)[&type=0]:0.5185631138)[&type=0]:0.3931733095)[&type=0]:0.4522656605)[&type=0]:0.5629860586):0.4667462698;";
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

        lifetime = new RealParameter("0.5 1 2 3");
        shape = new IntegerParameter("100 50 20 5");
        death = new RealParameter("0.05 0.1 0.1 0.1");
        sampling = new RealParameter("0.1 0.1 0.1 0.1");

        model = new Parameterization();
        model.initByName("nTypes", 4,
                "lifetime", lifetime,
                "shape", shape,
                "death", death,
                "symTransitions", new RealParameter("0.4 0.1 0 0 0 0.2 0 0 0 0 1 0 0 0 0 1"),
                "asymTransitions", new RealParameter("0 0.5 0 0 0 0 0.5 0.3 0 0 0 0 0 0 0 0"),
                "sampling", sampling,
                "originTime", 5.0);
    }

    @AfterEach
    public void printResult() {
        System.out.println(logHR);
        if (logHR != Double.NEGATIVE_INFINITY) {
            flatTree = tree.convertAnnotatedTree(true);
            System.out.println(flatTree.getRoot().toNewick());
        }
    }

    @Test
    public void testEventSampler() {

        // initialize operator
        EventSampler operator = new EventSampler();
        //lifetime.setValue(1, 0.5);
        operator.initByName("tree", tree, "parameterization", model, "drawEventCount", true, "proportionBranches", 1.0, "weight", 1.0);

        // trigger proposal
        logHR = operator.proposal();
    }


    @Test
    public void testWilsonBalding() {

        // initialize tree
        String newick = "(((((((((((((((((2[&type=2]:1.028989621)[&type=1]:1.013546413)[&type=1]:0.9258614777)[&type=0]:0.3843474571)[&type=0]:0.5182004427)[&type=0]:0.5348274266)[&type=0]:0.5238443404)[&type=0]:0.4389043469)[&type=0]:0.479401034)[&type=0]:0.4575102818)[&type=0]:0.4532865737)[&type=0]:0.4466242398)[&type=0]:0.4251882783)[&type=0]:0.4977463013)[&type=0]:0.4500339293,((((((((((((((((5[&type=1]:0.4468829172)[&type=0]:0.4496776294)[&type=0]:0.4669761298)[&type=0]:0.4961852113)[&type=0]:0.6053564319)[&type=0]:0.5116818928)[&type=0]:0.5593359388)[&type=0]:0.5041445555)[&type=0]:0.4352931526)[&type=0]:0.4451114997)[&type=0]:0.5905145662)[&type=0]:0.5563581743)[&type=0]:0.4359201852)[&type=0]:0.4723996884)[&type=0]:0.4850848771,((((((((((((4[&type=1]:0.5807148434)[&type=1]:1.057636084)[&type=0]:0.5225289534)[&type=0]:0.5476839539)[&type=0]:0.4806513184)[&type=0]:0.5873693377)[&type=0]:0.4439127757)[&type=0]:0.5123293585)[&type=0]:0.5599312644)[&type=0]:0.6005180287)[&type=0]:0.485795349)[&type=0]:0.4848830855,((((((((1[&type=2]:1.618415343)[&type=2]:1.32887674,((3[&type=3]:1.024651724)[&type=1]:0.9009859527)[&type=1]:1.021654405)[&type=1]:0.8688355805)[&type=0]:0.6204535167)[&type=0]:0.4382648793)[&type=0]:0.5024736266)[&type=0]:0.4562760577)[&type=0]:0.5162409957)[&type=0]:0.5141176133)[&type=0]:0.596968498)[&type=0]:0.618101857)[&type=0]:0.4992874571)[&type=0]:0.4856169175,(((((6[&type=2]:0.1827961005)[&type=2]:1.918497782)[&type=2]:2.430571703)[&type=2]:2.125177571)[&type=2]:1.519962828)[&type=1]:0.8869230981)[&type=0]:0.5037850172):0.4537720971;" ;
        AnnotatedTreeParser tree = new AnnotatedTreeParser();
        tree.initByName("newick", newick);

        // initialize operator
        AnnotatedWilsonBalding operator = new AnnotatedWilsonBalding();
        operator.initByName("tree", tree, "weight", 1.0);

        // trigger proposal
        operator.proposal();
    }


    @Test
    public void testExchange() {

        // initialize tree
        String newick = "(((((((((((((((((2[&type=2]:1.028989621)[&type=1]:1.013546413)[&type=1]:0.9258614777)[&type=0]:0.3843474571)[&type=0]:0.5182004427)[&type=0]:0.5348274266)[&type=0]:0.5238443404)[&type=0]:0.4389043469)[&type=0]:0.479401034)[&type=0]:0.4575102818)[&type=0]:0.4532865737)[&type=0]:0.4466242398)[&type=0]:0.4251882783)[&type=0]:0.4977463013)[&type=0]:0.4500339293,((((((((((((((((5[&type=1]:0.4468829172)[&type=0]:0.4496776294)[&type=0]:0.4669761298)[&type=0]:0.4961852113)[&type=0]:0.6053564319)[&type=0]:0.5116818928)[&type=0]:0.5593359388)[&type=0]:0.5041445555)[&type=0]:0.4352931526)[&type=0]:0.4451114997)[&type=0]:0.5905145662)[&type=0]:0.5563581743)[&type=0]:0.4359201852)[&type=0]:0.4723996884)[&type=0]:0.4850848771,((((((((((((4[&type=1]:0.5807148434)[&type=1]:1.057636084)[&type=0]:0.5225289534)[&type=0]:0.5476839539)[&type=0]:0.4806513184)[&type=0]:0.5873693377)[&type=0]:0.4439127757)[&type=0]:0.5123293585)[&type=0]:0.5599312644)[&type=0]:0.6005180287)[&type=0]:0.485795349)[&type=0]:0.4848830855,((((((((1[&type=2]:1.618415343)[&type=2]:1.32887674,((3[&type=3]:1.024651724)[&type=1]:0.9009859527)[&type=1]:1.021654405)[&type=1]:0.8688355805)[&type=0]:0.6204535167)[&type=0]:0.4382648793)[&type=0]:0.5024736266)[&type=0]:0.4562760577)[&type=0]:0.5162409957)[&type=0]:0.5141176133)[&type=0]:0.596968498)[&type=0]:0.618101857)[&type=0]:0.4992874571)[&type=0]:0.4856169175,(((((6[&type=2]:0.1827961005)[&type=2]:1.918497782)[&type=2]:2.430571703)[&type=2]:2.125177571)[&type=2]:1.519962828)[&type=1]:0.8869230981)[&type=0]:0.5037850172):0.4537720971;" ;
        AnnotatedTreeParser tree = new AnnotatedTreeParser();
        tree.initByName("newick", newick);

        // initialize operator
        AnnotatedExchange operator = new AnnotatedExchange();
        operator.initByName("tree", tree, "weight", 1.0, "isNarrow", false);

        // trigger proposal
        operator.proposal();
        Tree flatTree = ((AnnotatedTree)tree).convertAnnotatedTree(true);
        System.out.println(flatTree.getRoot().toNewick());
    }


    @Test
    public void testSubtreeSlide() {

        // initialize tree
        String newick = "(((((((((((((((((2[&type=2]:1.028989621)[&type=1]:1.013546413)[&type=1]:0.9258614777)[&type=0]:0.3843474571)[&type=0]:0.5182004427)[&type=0]:0.5348274266)[&type=0]:0.5238443404)[&type=0]:0.4389043469)[&type=0]:0.479401034)[&type=0]:0.4575102818)[&type=0]:0.4532865737)[&type=0]:0.4466242398)[&type=0]:0.4251882783)[&type=0]:0.4977463013)[&type=0]:0.4500339293,((((((((((((((((5[&type=1]:0.4468829172)[&type=0]:0.4496776294)[&type=0]:0.4669761298)[&type=0]:0.4961852113)[&type=0]:0.6053564319)[&type=0]:0.5116818928)[&type=0]:0.5593359388)[&type=0]:0.5041445555)[&type=0]:0.4352931526)[&type=0]:0.4451114997)[&type=0]:0.5905145662)[&type=0]:0.5563581743)[&type=0]:0.4359201852)[&type=0]:0.4723996884)[&type=0]:0.4850848771,((((((((((((4[&type=1]:0.5807148434)[&type=1]:1.057636084)[&type=0]:0.5225289534)[&type=0]:0.5476839539)[&type=0]:0.4806513184)[&type=0]:0.5873693377)[&type=0]:0.4439127757)[&type=0]:0.5123293585)[&type=0]:0.5599312644)[&type=0]:0.6005180287)[&type=0]:0.485795349)[&type=0]:0.4848830855,((((((((1[&type=2]:1.618415343)[&type=2]:1.32887674,((3[&type=3]:1.024651724)[&type=1]:0.9009859527)[&type=1]:1.021654405)[&type=1]:0.8688355805)[&type=0]:0.6204535167)[&type=0]:0.4382648793)[&type=0]:0.5024736266)[&type=0]:0.4562760577)[&type=0]:0.5162409957)[&type=0]:0.5141176133)[&type=0]:0.596968498)[&type=0]:0.618101857)[&type=0]:0.4992874571)[&type=0]:0.4856169175,(((((6[&type=2]:0.1827961005)[&type=2]:1.918497782)[&type=2]:2.430571703)[&type=2]:2.125177571)[&type=2]:1.519962828)[&type=1]:0.8869230981)[&type=0]:0.5037850172):0.4537720971;";
        AnnotatedTreeParser tree = new AnnotatedTreeParser();
        tree.initByName("newick", newick);

        // initialize operator
        AnnotatedSubtreeSlide operator = new AnnotatedSubtreeSlide();
        operator.initByName("tree", tree, "weight", 1.0, "size", 3);

        // trigger proposal
        operator.proposal();
    }

}
