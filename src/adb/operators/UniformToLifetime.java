/*
* File Uniform.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
/*
 * UniformOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package adbp.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;


@Description("Randomly selects true internal tree node (i.e. not the root) and move node height uniformly in interval " +
        "restricted by the nodes parent and children.")
public class UniformToLifetime extends TreeOperator {

    final public Input<RealParameter> expectedLifetime = new Input<>("expectedLifetime", "expected lifetime");
    final public Input<Double> error = new Input<>("error", "expected lifetime error");
    final public Input<Double> proportionInput = new Input<>("prop", "what is maximum percentage of suitable branches to rescale");

    double expectedLifetimeValue;
    double errorValue;
    double minLength = 0.0;
    private double logHR;

    @Override
    public void initAndValidate() {
        expectedLifetimeValue = expectedLifetime.get().getValue();
        errorValue = error.get();
        minLength = expectedLifetimeValue + errorValue;
        
        // check if the expected lifetime is smaller than the error
        if (expectedLifetimeValue < errorValue) {
            throw new IllegalArgumentException("Expected lifetime (" + expectedLifetimeValue + ") should be larger than error (" + errorValue + ")");
        }
        logHR = 0;
        
    }

    /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final Tree tree = (Tree) InputUtil.get(treeInput, this);
        logHR = 0;

        
        // Abort if no non-root internal nodes
        if (tree.getInternalNodeCount()==1)
            return Double.NEGATIVE_INFINITY;

        List<Node> canScaleBefore = canScale(minLength);
        
        if (canScaleBefore.isEmpty()) {
            return Double.NEGATIVE_INFINITY;
        }
        
        // randomly select internal node
        List<Node> nodes = new ArrayList<>();
        do {
            Node node;
            do{
                node = canScaleBefore.get(Randomizer.nextInt(canScaleBefore.size()));
            } while (node.isRoot() || node.isLeaf());
            nodes.add(node);
        } while ((nodes.size()==1 && proportionInput.get()==null) ||
                (double) nodes.size() / canScaleBefore.size() < proportionInput.get());


        Collections.sort(nodes, Comparator.comparing(u -> u.getHeight()));

        if (Randomizer.nextDouble() < 0.5) {
            for (int i = nodes.size(); i>0; i--){
                Node node = nodes.get(i-1);

                double newValue;
                double upper = node.getParent().getHeight();
                double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());

                // move to the mean
                // double upper_ = upper - expectedLifetimeValue + errorValue;
                double lower_ = upper - expectedLifetimeValue - errorValue;
                newValue =  (Randomizer.nextDouble() * (2*errorValue)) + lower_ ;
                if (newValue < lower || newValue > upper) {
                    // invalid move, can be rejected immediately
                    return Double.NEGATIVE_INFINITY;
                }

                // this action
                logHR -= Math.log(1/(2*errorValue));
                // and the reverse action
                logHR += Math.log(1/(upper - lower));
                node.setHeight(newValue);
            }
        } else {
            for (Node node : nodes){
                double newValue;
                double upper = node.getParent().getHeight();
                double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
                newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;
                if (newValue < lower || newValue > upper) {
                    // invalid move, can be rejected immediately
                    return Double.NEGATIVE_INFINITY;
                }
                // reverse action
                logHR += Math.log(1/(2*errorValue));
                // this action
                logHR -= Math.log(1/(upper - lower));
                node.setHeight(newValue);

            }
        }


        List<Node> canScaleAfter = canScale(minLength);
        logHR += Math.log(nodes.size()/canScaleAfter.size());
        logHR -= Math.log(nodes.size()/canScaleBefore.size());

        
        

        return logHR;
    }

    private List<Node> canScale(double minLength){
        List<Node> nodes = new ArrayList<>();
        for (Node node : treeInput.get().getInternalNodes()) {
            if (node.isRoot()) {
                continue;
            }
            final double upper = node.getParent().getHeight();
            final double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
            if ((upper - lower) >= minLength) {
                nodes.add(node);
            }
        }
        return nodes;
    }

}
