package test.adbp;

import multitypetree.evolution.tree.MultiTypeNode;
import multitypetree.evolution.tree.MultiTypeTreeFromNewick;
import org.junit.jupiter.api.Test;

public class MultiTypeTreeTest {

    @Test
    public void testMultiTypeTree() {

        String newick = "(1[&state=0] : 1.5, 2[&state=1] : 0.5)[&state=0];";

        MultiTypeTreeFromNewick mtt = new MultiTypeTreeFromNewick();
        mtt.initByName(
                "adjustTipHeights", false,
                "value", newick,
                "typeLabel", "state");

        int n = mtt.getTypeSet().getNTypes();
        int m = mtt.getNodeCount();

        for (int i = 0; i < m; i++) {
            MultiTypeNode node = (MultiTypeNode) mtt.getNode(i);
            System.out.println(node.getNodeType());
        }
    }
}
