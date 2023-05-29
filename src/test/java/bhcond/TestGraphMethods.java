package bhcond;

import bhcond.BasicDendrogram;
import bhcond.BasicTreeNode;
import bhcond.GraphAdjMatrix;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

import com.google.common.collect.Iterables;
import com.google.common.graph.*;

public class TestGraphMethods {
	
	public static void main(String[] args) {
		System.out.println("~~Test Program for Dendrogram and Graph Methods~~ \n");
		
		// Create the nodes
		BasicTreeNode root = BasicTreeNode.createInternal(1, null);
		BasicTreeNode leafLeft = BasicTreeNode.createLeaf(1, true);
		BasicTreeNode leafRight = BasicTreeNode.createLeaf(2, false);
		
		// Build a small dendrogram with 2 leaves
		BasicDendrogram dend = buildSmallDend(root, leafLeft, leafRight);
	
		// Run some checks on the core dendrogram methods
		runBasicChecks(dend, root, leafLeft, leafRight);
		
		///////////////////////////////////
		// Experiment on a larger example:
		///////////////////////////////////
		
		// Build a dendrogram with n leaves
		int numLeaves = 4;
		int[] vertices = IntStream.rangeClosed(1, numLeaves).toArray();
		//System.out.println(Arrays.toString(vertices));
		GraphAdjMatrix graph = new GraphAdjMatrix(vertices);
		BasicDendrogram d = BasicDendrogram.initDend(graph);
		// note: recall that initDend performs a shuffle of vertices
		
		// Establish the structure of this new dendrogram and experiment with interchange moves
		
		System.out.println(d);
		System.out.println(d.prettyPrint());
		
		BasicTreeNode top = d.getRoot();
		List<Integer> leaves = d.getLeaves(top);
		
		//System.out.println(d.printLeavesPartition(top));
		
		System.out.println("the leaves are: " + leaves);
		
		
		System.out.println();
		
		// get the internal nodes
		Set<BasicTreeNode> internalNodes = d.getInternalNodes(); // includes root
		
		// remove the root
		Iterator<BasicTreeNode> it = internalNodes.iterator();
		BasicTreeNode nextNode = it.next();
		while(it.hasNext() && nextNode.isLeftChild() != null) {
			nextNode = it.next();
		}
		internalNodes.remove(nextNode);
		
		// choose an internal node
		it = internalNodes.iterator();
		BasicTreeNode internal = it.next();
		System.out.println("The chosen internal node is: " + internal);
		System.out.println("The chosen internal node is a left child: " + internal.isLeftChild());
		
		// choose a child node
		Set<BasicTreeNode> children = d.getChildren(internal);
		Iterator<BasicTreeNode> iter = children.iterator();
		BasicTreeNode child = iter.next();
		System.out.println("The chosen child node is: " + child + " with vertex ID " + child.getVertexID());
		System.out.println("The chosen child node is a left child: " + child.isLeftChild());
		
		// perform interchange move
		System.out.println("Performing an interchange move...");
		d.interchange(internal, child);
		
		System.out.println(d.prettyPrint());
		System.out.println("The dendrogram viewed as a partition of the leaves:");
		System.out.println(d);
		System.out.println();
		System.out.println("The chosen internal node is a left child: " + internal.isLeftChild());
		System.out.println("The chosen child node is a left child: " + child.isLeftChild());
		
		// revert the interchange move and inspect
		System.out.println("Reverting the interchange...");
		d.interchange(internal, child);
		
		System.out.println(d.prettyPrint());
		System.out.println("The dendrogram viewed as a partition of the leaves:");
		System.out.println(d);
		System.out.println();
		System.out.println("The chosen internal node is a left child: " + internal.isLeftChild());
		System.out.println("The chosen child node is a left child: " + child.isLeftChild());

		
		
	}
	
	public static void runBasicChecks(BasicDendrogram dend, BasicTreeNode root, BasicTreeNode leafLeft, BasicTreeNode leafRight) {
		
		System.out.println(dend);
		//System.out.println("The dendrogram viewed as a partition of the leaves:");
		//System.out.println(dend.printLeavesPartition(root));
		System.out.println(dend.prettyPrint());
		
		// count leaves
		System.out.println("the number of leaves in the dendrogram is: " + dend.getNumLeaves());
		System.out.println("the number of leaves in the left subtree is: " + dend.countLRleaves(root, true));
		System.out.println("the number of leaves in the right subtree is: " + dend.countLRleaves(root, false) + "\n");
		
		// check if leaf or internal
		System.out.println("leafLeft is a leaf: " + dend.isLeaf(leafLeft));
		System.out.println("leafRight is a leaf: " + dend.isLeaf(leafRight));
		System.out.println("root is a NOT a leaf: " + !dend.isLeaf(root) + "\n");
		
		// check reachability 
		System.out.println("leafLeft is reachable from root: " + dend.isReachable(root, leafLeft));
		System.out.println("leafRight is reachable from root: " + dend.isReachable(root, leafRight));
		System.out.println("leafLeft is NOT reachable from leafRight: " + !dend.isReachable(leafRight, leafLeft));
		System.out.println("leafRight is NOT reachable from leafLeft: " + !dend.isReachable(leafLeft, leafRight));
		System.out.println("root is NOT reachable from leafLeft: " + !dend.isReachable(leafLeft, root));
		System.out.println("root is NOT reachable from leafRight: " + !dend.isReachable(leafRight, root) + "\n");
		
		// check vertex IDs and getChildren (successors)
		System.out.println("the vertex ID of leafLeft is: " + leafLeft.getVertexID());
		System.out.println("the vertex ID of leafRight is: " + leafRight.getVertexID());
		
		Set<BasicTreeNode> children = dend.getChildren(root);
		System.out.println("the vertex IDs of the children of root are: ");
		for (BasicTreeNode n : children) {
		    System.out.println(n.getVertexID());
		}
		System.out.println();
		
		// check getParent (predecessors)
		System.out.println("the parent of leafLeft has a vertex ID of: " +
				dend.getParent(leafLeft).getVertexID());
		System.out.println("the parent of leafLeft is NOT a leaf: " +
				!dend.isLeaf(dend.getParent(leafLeft)));
		System.out.println("the parent of leafLeft has a left property of: " +
				dend.getParent(leafLeft).isLeftChild() + "\n");
		
		System.out.println("the parent of leafRight has a vertex ID of: " +
				dend.getParent(leafRight).getVertexID());
		System.out.println("the parent of leafRight is NOT a leaf: " +
				!dend.isLeaf(dend.getParent(leafRight)));
		System.out.println("the parent of leafRight has a left property of: " +
				dend.getParent(leafRight).isLeftChild() + "\n");
		
		// check getSubtree (reachableNodes)
		Set<BasicTreeNode> subtree = dend.getSubtree(root);
		System.out.println("the vertex IDs of the subtree rooted from 'root' are: ");
		for (BasicTreeNode n : subtree) {
		    System.out.println(n.getVertexID());
		}
		
		// check getLRSubtree
		List<Integer> leftSubtree = dend.getLRSubtree(root, true);
		System.out.println("the vertex ID(s) of the leaves from the left subtree rooted from 'root' are: " + 
				Arrays.toString(leftSubtree.toArray()));
		System.out.println("the number of leaves on the left subtree rooted from 'root' are: " + 
				dend.countLRleaves(root, true));
		List<Integer> rightSubtree = dend.getLRSubtree(root, false);
		System.out.println("the vertex ID(s) of the leaves from the right subtree rooted from 'root' are: " + 
				Arrays.toString(rightSubtree.toArray()));
		System.out.println("the number of leaves on the right subtree rooted from 'root' are: " + 
				dend.countLRleaves(root, false));
		List<Integer> leftSubtree2 = dend.getLRSubtree(leafLeft, true);
		System.out.println("the vertex ID(s) of the leaves from the left subtree rooted from 'leafLeft' are: " + 
				Arrays.toString(leftSubtree2.toArray()) + "\n");
		
		// check getInternalNodes
		Set<BasicTreeNode> internal = dend.getInternalNodes();
		System.out.println("the vertex IDs of the internal nodes are: ");
		for (BasicTreeNode n : internal) {
		    System.out.println(n.getVertexID());
		}
		
		// ************************
		// check interchange method - requires a tree with at least 3 leaves
		//dend.interchange(root, leafLeft); // should throw Runtime Exception
		
		System.out.println();
	}
	
	public static BasicDendrogram buildSmallDend(BasicTreeNode root, BasicTreeNode leafLeft, BasicTreeNode leafRight) {
				
		// Create the dendrogram
		BasicDendrogram dend = new BasicDendrogram();
				
		// Construct the dendrogram
		dend.setRoot(root);
		dend.addEdge(root, leafLeft);
		dend.addEdge(root, leafRight);
		
		return(dend);
	}
	

}
