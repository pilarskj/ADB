package adb.operators;

import adb.tree.AnnotatedTree;
import adb.tree.AnnotatedTreeParser;
import beast.base.evolution.tree.Tree;
import org.junit.jupiter.api.Test;

public class MTAnnotatedTreeOperatorTest {

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
