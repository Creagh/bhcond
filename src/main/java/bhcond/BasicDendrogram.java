/**
 * 
 */
package bhcond;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

import blang.inits.experiments.tabwriters.TidilySerializable;
import blang.inits.experiments.tabwriters.TidySerializer.Context;
import blang.mcmc.Samplers;
import blang.runtime.internals.objectgraph.DeepCloner;
import briefj.collections.UnorderedPair;

import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.graph.*;
import com.google.common.primitives.Ints;
import com.rits.cloning.IDeepCloner;
import com.rits.cloning.IFastCloner;

import bayonet.math.SamplingUtils;
import bhcond.BasicTreeNode;

/**
 * A dendrogram is a rooted binary tree. We represent a tree as a directed graph.
 * 
 * FOR USE WITH GraphAdjMatrix and BasicTreeNode classes
 * 
 * The annotation "@Samplers" links the data type with the appropriate sampler. 
 * 
 * 
 * @author Creagh Briercliffe (creagh.briercliffe@gmail.com)
 *
 */
@Samplers(DendrogramSampler.class)
public class BasicDendrogram {
  
  /*
   * For unknown reason, automatic cloning of Google Graph objects via DeepCloner fails.
   * This is a workaround which manually create deep copies.
   */
  static {
    DeepCloner.cloner.registerFastCloner(BasicDendrogram.class, new IFastCloner() {
      @SuppressWarnings({ "rawtypes", "unchecked" })
      @Override
      public Object clone(Object t, IDeepCloner cloner, Map<Object, Object> clones) {
        BasicDendrogram input = (BasicDendrogram) t;
        MutableGraph copy = GraphBuilder.directed().build();
        for (Object e : input.tree.edges()) {
          // deep clone the edge
          EndpointPair eCopy = (EndpointPair) cloner.deepClone(e, clones);
          copy.putEdge(eCopy.source(), eCopy.target());
        }
        BasicTreeNode rootCopy = (BasicTreeNode) cloner.deepClone(input.root, clones);
        return new BasicDendrogram(copy, rootCopy);
      }
    });
  }

	private MutableGraph<BasicTreeNode> tree;
	private BasicTreeNode root; // the root node of the dendrogram

	public BasicDendrogram(MutableGraph<BasicTreeNode> tree, BasicTreeNode root) {
		this.tree = tree;
		this.root = root;
	}

	/**
	 * Create an empty dendrogram as a directed graph
	 */
	public BasicDendrogram() {
		tree = GraphBuilder.directed().build();
	}
	
	/**
	 * Set the tree portion of the dendrogram
	 */
	void setTree(MutableGraph<BasicTreeNode> tree) {
		this.tree = tree;
	}
	
	/**
	 * Set the root node of the tree
	 * 
	 * @param root
	 */
	void setRoot(BasicTreeNode top) {
		if (top == null)
			throw new RuntimeException();
		this.root = top;
	}
	
	/**
	 * Get the root of the dendrogram
	 * 
	 */
	BasicTreeNode getRoot() {
		return this.root;
	}
	
	/**
	 * Get the tree portion of the dendrogram
	 * 
	 */
	MutableGraph<BasicTreeNode> getTree(){
		return this.tree;
	}

	/**
	 * Add a node to the tree
	 * 
	 * @param node
	 */
	void addNode(BasicTreeNode node) {
		this.tree.addNode(node);
	}

	/**
	 * Add a directed edge from node 1 to node 2, if one is not already present
	 * 
	 * @param node1
	 * @param node2
	 */
	void addEdge(BasicTreeNode node1, BasicTreeNode node2) {
		
		// check to ensure that there isn't already an edge from node 2 to node 1
		//if (tree.hasEdgeConnecting(node2, node1)) {
		//    throw new IllegalArgumentException("Cannot add this edge because node2 already has an edge to node1.");  
		//} else {
			this.tree.putEdge(node1, node2);
		//}
	}
	
	/**
	 * Check to see if a target node is reachable from a given node
	 * 
	 * @param from
	 * @param target
	 * @return boolean true if node from can reach node target
	 */
	boolean isReachable(BasicTreeNode from, BasicTreeNode target) {
		
		boolean canReach = false;
		Set<BasicTreeNode> reachables = Graphs.reachableNodes(this.tree, from);
		
		// iterate through the reachables set to see if it contains target node
		Iterator<BasicTreeNode> it = reachables.iterator();
		while(!canReach && it.hasNext()) {
			BasicTreeNode nextNode = it.next();
			if(nextNode == target) { // target is reachable
				canReach = true;
			}
		}
		
		return canReach;
	}

	/**
	 * Remove the edge from node 1 to node 2, if it is present
	 * 
	 * @param node1
	 * @param node2
	 */
	void removeEdge(BasicTreeNode node1, BasicTreeNode node2) {
		this.tree.removeEdge(node1, node2);
	}
	
	/**
	 * Get the parent of a node; return null if node is root.
	 * 
	 * @param node
	 * @return parent
	 */
	public BasicTreeNode getParent(BasicTreeNode node) {
		Set<BasicTreeNode> predecessors = this.tree.predecessors(node);
		BasicTreeNode parent = null;
		
		if(predecessors.size() == 1) { // a non-root node in a binary tree must have only a single parent
			//parent = predecessors.iterator().next();
			parent = Iterables.getOnlyElement(predecessors);

		} else if(predecessors.size() > 1) { 
			throw new RuntimeException("A node of a binary tree should not have more than one parent.");
		}
		return(parent);
	}
	
	/**
	 * Get the children of a node (does not include grandchildren, etc.)
	 * 
	 * @param node
	 * @return children
	 */
	public Set<BasicTreeNode> getChildren(BasicTreeNode node) {
		return(this.tree.successors(node));
	}
	
	/**
	 * Get the subtree rooted at a given node. This
	 * returns all descendants of a node, including the original node.
	 * 
	 * @param node
	 * @return descendants
	 */
	public Set<BasicTreeNode> getSubtree(BasicTreeNode node) {
		return(Graphs.reachableNodes(this.tree, node));
	}
	
	/**
	 * Determine if a given node from this dendrogram is a leaf
	 * 
	 * @param node
	 * @return boolean leaf
	 */
	public boolean isLeaf(BasicTreeNode node) {
		boolean leaf = false;
		
		if(this.tree.outDegree(node) == 0) { // zero outgoing edges
			leaf = true;
		}
		return(leaf);
	}
	
	/**
	 * Get the vertex ID's of the leaves from the subtree rooted at a given node.
	 * If the given node is a leaf, then this function returns the given node's ID.
	 * 
	 * @param node
	 * @return leaves
	 */
	public List<Integer> getLeaves(BasicTreeNode node) {
		
		List<Integer> leaves = new ArrayList<Integer>();
		
		// TODO: determine which version of isLeaf is faster, which is also used in getLRSubtree() method
		
		// check to see if node is a leaf
		if(isLeaf(node)) {
		//if(node.isLeaf()) {
			leaves.add(node.getVertexID());
		} else { 
			// move down the tree to find the leaves:
			// get the subtree rooted at node (which includes node)
			Set<BasicTreeNode> subtree = getSubtree(node);
			// construct an iterator over the elements of subtree
			Iterator<BasicTreeNode> it = subtree.iterator();
			
			// NOTE: subtree includes 'node', which we already know isn't a leaf;
			// hence, we do one more check than necessary. Instead we could use subtree.remove(node),
			// but that would require a method to check equality of nodes.

			// iterate through the subtree to find the leaves
			while(it.hasNext()) {
				BasicTreeNode nextNode = it.next();
				if(isLeaf(nextNode)) { // nextNode is a leaf
				//if(nextNode.isLeaf()) {
					leaves.add(nextNode.getVertexID());
				}
			}
		}
		
		return(leaves);
	}
	
	/**
	 * Get the number of leaves in the tree
	 * 
	 * @return int numLeaves
	 */
	public int getNumLeaves() {
		List<Integer> leaves = getLeaves(this.root);
		
		return(leaves.size());
	}
	
	/**
	 * Get the vertex ID's of the leaves from the specified left or right subtree rooted at a given node.
	 * If the given node is a leaf, then this function returns the given node.
	 * 
	 * @param node
	 * @param left - boolean TRUE if left subtree; FALSE if right subtree
	 * @return leaves - list of integer ID's for leaves
	 */
	public List<Integer> getLRSubtree(BasicTreeNode node, boolean left) {
		
		List<Integer> leaves = new ArrayList<Integer>();
		
		// check to see if node is a leaf
		if(isLeaf(node)) {
			leaves.add(node.getVertexID());
		} else { 
			// Find the appropriate subtree (left/right):
			// first, get the children of node
			Set<BasicTreeNode> children = getChildren(node);
			
			// construct an iterator over the elements of children
			Iterator<BasicTreeNode> it = children.iterator();
			// iterate through the children to find the appropriate (left/right) child
			BasicTreeNode nextNode = it.next();
			if(nextNode.isLeftChild() == null) {
				//System.out.println("nextNode is root...");
				throw new RuntimeException("A child node should not have a left property that is null.");
			}
			if(left != nextNode.isLeftChild()) { // nextNode is not the appropriate subtree (left/right)
				nextNode = it.next(); // iterate to the only other (and therefore appropriate) child
			} 
			
			// get the leaves from the subtree rooted at nextNode
			leaves = getLeaves(nextNode);	

		}
		
		return leaves;
	}
	
	
	/**
	 * Get the set of all internal nodes in the dendrogram
	 * i.e. all nodes (including root) except for the leaves
	 * 
	 * @return a set of a BasicTreeNode objects
	 */
	public Set<BasicTreeNode> getInternalNodes() {
		BasicTreeNode root = this.root;
		
		// convert the tree to a set of nodes (as opposed to a mutable graph)
		// by extracting the subtree rooted at the root node
		Set<BasicTreeNode> treeSet = new LinkedHashSet<BasicTreeNode>(getSubtree(root));
		
		// iterate through the treeSet to find and remove the leaves
		for (Iterator<BasicTreeNode> it = treeSet.iterator(); it.hasNext();) {
			BasicTreeNode nextNode = it.next();
		    if (isLeaf(nextNode)) { // nextNode is a leaf
		        it.remove(); // pluck the leaf
		    }
		}
		
		return(treeSet);
		
	}
	
	/**
	 * Count the number of leaves from the specified left or right subtree rooted at a given node.
	 * If the given node is a leaf, then this function returns 1.
	 * 
	 * (i.e. count L_r or R_r for an internal node r)
	 * 
	 * @param node
	 * @param left - boolean TRUE if left subtree; FALSE if right subtree
	 * @return numLeaves 
	 */
	public int countLRleaves(BasicTreeNode node, Boolean left) {
		List<Integer> leaves = getLRSubtree(node, left);
		return(leaves.size());
	}
	
	
	/**
	 * Perform a Nearest Neighbour Interchange type move for a rooted tree.
	 * 
	 * This method rearranges the subtrees in a dendrogram as follows.
	 * For an internal node, r, there are two subtrees, s and t, descended
	 * from r's children. There's a third subtree, u, descended from r's sibling.
	 * There are two ways to rearrange these subtrees without disturbing their internal
	 * structure. 
	 * 
	 * 
	 * These two rearrangements can be constructed with a call to either 
	 * interchange(r, s) or interchange(r, t).
	 * 
	 * @param node - an internal node (not the root node)
	 * @param child - one of node's children
	 */
	public void interchange(BasicTreeNode node, BasicTreeNode child) {
		// Let node = r and child = s.
		// Let child's sibling be denoted by t and child's uncle be denoted by u.
		// Let node's parent be denoted by q.
		
		// the result of calling interchange(r, s)
		//     q                q 
		//    / \              / \
		//   r   \    ===>    /   r				(1)
		//  / \   \          /   / \
		// s   t   u        t   s   u
		
		//     q                q 
		//    / \              / \
		//   r   \    ===>    /   r				(2)
		//  / \   \          /   / \
		// t   s   u        t   u   s
		
		//     q                q 
		//    / \              / \
		//   /   r    ===>    r   \				(3)
		//  /   / \          / \   \
		// u   s   t        s   u   t
		
		//     q                q 
		//    / \              / \
		//   /   r    ===>    r   \				(4)
		//  /   / \          / \   \
		// u   t   s        u   s   t
		
		// The LHS of (1) - (4) all have the same topology 
		// (and are only distinct when we keep track of left and right properties).
		// Similarly the RHS of (1) - (4) all share the same topology.

		// the result of calling interchange(r, t)
		//     q                q 
		//    / \              / \
		//   r   \    ===>    /   r				(5)
		//  / \   \          /   / \
		// s   t   u        s   u   t
		
		//     q                q 
		//    / \              / \
		//   r   \    ===>    /   r				(6)
		//  / \   \          /   / \
		// t   s   u        s   t   u
		
		//     q                q 
		//    / \              / \
		//   /   r    ===>    r   \				(7)
		//  /   / \          / \   \
		// u   s   t        u   t   s
		
		//     q                q 
		//    / \              / \
		//   /   r    ===>    r   \				(8)
		//  /   / \          / \   \
		// u   t   s        t   u   s
		
		// The RHS of (5) - (8) all share the same topology, which is distinct from the RHS of (1) - (4).
		
		// Given that BasicTreeNode records a left/right property, we take the following convention.
		
		// CONVENTION:
		// interchange(node, child) swaps the left/right property of node;
		// child retains their left/right property;
		// if node's new property matches former sibling of child, then swaps sibling's left/right property;
		// former uncle takes the opposite property of child;
		
		//if(getParent(node) == null) {
		if(node.isLeftChild() == null) {
			throw new RuntimeException("Cannot perform an interchange for the root node.");
		} else {
			
			// get the left/right property of node
			boolean nodeIsLeft = node.isLeftChild();
			
			// swap the left/right property of node
			nodeIsLeft = !nodeIsLeft;
			// set the left/right property of node as the swapped value
			node.setLeft(nodeIsLeft);
			
			// remove the edge between node and sibling of child:
			// this is done by temporarily removing the edge from node to child
			removeEdge(node, child);
			// get the sibling of child
			BasicTreeNode sibling = getChildren(node).iterator().next();
			// remove the edge between node and sibling of child
			removeEdge(node, sibling);
			// add back the edge between node and child
			addEdge(node, child);

			// remove the edge between parent and uncle:
			// this is done by temporarily removing the edge from parent to node
			// in order to get uncle
			BasicTreeNode parent = getParent(node);
			removeEdge(parent, node);
			BasicTreeNode uncle = getChildren(parent).iterator().next();
			// remove the edge between parent and uncle
			removeEdge(parent, uncle);
			// add back the edge between parent and node
			addEdge(parent, node);

			// add an edge from parent to sibling of child
			addEdge(parent, sibling);

			// add an edge from node to uncle
			addEdge(node, uncle);
			
			// perform the appropriate left/right property switching according to the
			// convention stated above:
			
			// if node's new property matches sibling's, then swaps sibling's left/right property;
			boolean siblingIsLeft = sibling.isLeftChild();
			
			if(nodeIsLeft == siblingIsLeft) { // they match
				siblingIsLeft = !siblingIsLeft; // swap left/right property
				sibling.setLeft(siblingIsLeft);
			}
			
			// uncle takes the opposite property of child
			uncle.setLeft(!child.isLeftChild());
		}

	}
	
	
	/**
	 * Create an initial Dendrogram, based on the graph, used for initializing the sampler.
	 * 
	 * A basic dendrogram is created by shuffling the vertices, then building a complete(*) binary tree.
	 * (*)The tree is filled left to right, by making everything a left child by default.
	 * 
	 * 
	 * @param graph
	 * @return dend
	 */
	public static  BasicDendrogram initDend(GraphAdjMatrix graph) {
		
		BasicDendrogram dend = new BasicDendrogram();
		int[] vertices = graph.getVertices();
		
		// store the tree nodes in a queue
		Queue<BasicTreeNode> queue = new LinkedList<BasicTreeNode>(); 
		
		// Build the tree from the bottom up:
		// first, shuffle the vertices
		shuffleArray(vertices);
		
		// create a leaf node for each vertex, and add the nodes to the queue
		for(int v : vertices) {

			Boolean left = true; // initially set all leaf nodes as left children
			
			// create the leaf node
			BasicTreeNode leaf = BasicTreeNode.createLeaf(v, left);
			// add the leaf node to the 	queue
			queue.add(leaf);
		}
		
		// while there is at least 2 nodes in the queue,
		// take a pair of nodes from the queue and connect with a new internal node
		// add the internal node to the queue.
		
		int id = 1;
		while(queue.size() > 1) {
			
			BasicTreeNode node1 = queue.poll();
			BasicTreeNode node2 = queue.poll();
			
			// make node2 the right child
			node2.setLeft(false);
			
			// create a new internal node
			Boolean left = true; // initially set all internal nodes as left children
			BasicTreeNode internal = BasicTreeNode.createInternal(id++, left);

			
			// add a directed edge from the internal node to nodes 1 & 2
			// NOTE: from Google Guava "You can add an edge whose incident nodes 
			// have not previously been added to the graph. If they're not already present, 
			// they're silently added to the graph"
			dend.addEdge(internal, node1);
			dend.addEdge(internal, node2);
			
			// add the internal node back to the queue
			queue.add(internal);
		}
		
		// The last node in the queue is the root node
		// Change the left property for the root node
		BasicTreeNode root = queue.poll();
		root.setLeft(null);
		// Set the root node
		dend.setRoot(root);
		
		//System.out.println("The initial dendrogram: \n" + dend.prettyPrint());
		
		return dend;
	}
	
	/**
	 * Shuffle the elements of an integer array, in place.
	 * 
	 * This method uses the Durstenfeld Shuffle, a modern version
	 * of the Fisher-Yates shuffle algorithm.
	 * See https://en.wikipedia.org/wiki/Fisher-Yates_shuffle#The_modern_algorithm
	 * 
	 * @param array
	 */
	private static void shuffleArray(int[] array) {
		
		int index, temp;
		Random rand = new Random(1);
		
		for(int i = array.length - 1; i > 0; i--) {
			
			index = rand.nextInt(i + 1);
			// swap elements
			temp = array[index];
			array[index] = array[i];
			array[i] = temp;
		}
		
	}
	
	/**
	 * Print the leaves of a subtree, rooted from a given node, as a partition
	 * 
	 * @param node - the top of the subtree
	 */
	public String printLeavesPartition(BasicTreeNode node) {
		if(node.isLeaf()) {
			return "" + node.getVertexID();
		} else {
			// recurse through left and right subtrees
			Set<BasicTreeNode> children = getChildren(node);
			List<String> childStrings = new ArrayList<>();
			for(BasicTreeNode n: children) {
			  childStrings.add(printLeavesPartition(n));
			}
			return "(" + Joiner.on(",").join(childStrings) + ")";
		}

	}
	
	/**
	 * Print the dendrogram in the Newick tree format
	 */
    public String toString() {
      return printLeavesPartition(getRoot()) + ";";
    }
    
	/**
	 * Print the entire dendrogram to string
	 * 
	 * Reference: code adapted from Stack Overflow user "VasiliNovikov" answer
	 * https://stackoverflow.com/questions/4965335/how-to-print-binary-tree-diagram/8948691#8948691
	 */
    public String prettyPrint() {
      StringBuilder buffer = new StringBuilder(50);
      print(this.root, buffer, "", "");
      return buffer.toString();
    }

    private void print(BasicTreeNode node, StringBuilder buffer, String prefix, String childrenPrefix) {
    		String vertexID = String.valueOf(node.getVertexID());
    		//if(vertexID.equals("-1")) {
    		//	vertexID = node.toString();
    		//}
    		buffer.append(prefix);
        buffer.append(vertexID);
        buffer.append('\n');
        
        Set<BasicTreeNode> children = getChildren(node);
        
        for (Iterator<BasicTreeNode> it = children.iterator(); it.hasNext();) {
            BasicTreeNode next = it.next();
            if (it.hasNext()) {
                print(next, buffer, childrenPrefix + "├── ", childrenPrefix + "│   ");
            } else {
                print(next, buffer, childrenPrefix + "└── ", childrenPrefix + "    ");
            }
        }
    }
	
	  
	/**
	 * Sample a dendrogram independently and uniformly at random (IN PLACE)
	 * 
	 * Method: (1) sample a binary unrooted tree, with n+1 leaves, uniformly at random. Then (2)
	 * transform into a rooted binary tree with n leaves.
	 * 
	 * This approach uses the fact that there exists a bijective function mapping
	 * unrooted binary trees with n+1 leaves to the space of rooted binary trees with n
	 * leaves. Reference: "Phylogenetics" by Charles Semple & Mike Steel (Prop. 2.2.3)
	 */
	public void sampleUniform(Random rand) { 
		
		// Create a new dendrogram object, consisting of a tree
		BasicDendrogram dend = new BasicDendrogram();
		
		// Create a list to store the tree nodes
		List<BasicTreeNode> nodes = new ArrayList<BasicTreeNode>();
		// Start by extracting the leaves
		List<Integer> leaves = this.getLeaves(this.root);
		Collections.sort(leaves);
		
		// Create a BasicTreeNode for each leaf/vertex, and add the nodes to the list
		for(int leafID : leaves) {
			
			Boolean left = true; // initially set all leaf nodes as left children
			
			// create the leaf node
			BasicTreeNode leaf = BasicTreeNode.createLeaf(leafID, left);
			// add the leaf node to the list
			nodes.add(leaf);
		}
		
		// Find the maximum leafID in the data
		int max = Collections.max(leaves);
		
		// Add an extra leaf, which will get pruned during the transformation
		// back into a rooted tree
		BasicTreeNode extra = BasicTreeNode.createLeaf(max+1, true);
		nodes.add(extra);
		
		// Shuffle the tree nodes list
		Collections.shuffle(nodes, rand);
		
		////////////////////////////////////////////////////////////////////////////////
		// STAGE 1: Sample an unrooted binary tree uniformly at random from a shuffled
		// list of n+1 tree nodes. 
		////////////////////////////////////////////////////////////////////////////////
		
		// Approach: given an unrooted tree with k leaves, sample a tree with
		// k+1 leaves by the following iterative process. Pick an edge at random
		// and insert a new internal node along that edge connecting to a new leaf.
		
		// Create a queue of the leaf nodes
	    Queue<BasicTreeNode> queue = Lists.newLinkedList(nodes);
	    
	    // Base cases:
	    if (queue.isEmpty()) {return;} // no nodes
	    
	    BasicTreeNode leaf1 = queue.poll();
   		
	    if (queue.isEmpty()) {  // only 1 node
	    		dend.addNode(leaf1);
	    		return;
	    	}
	    
	    BasicTreeNode leaf2 = queue.poll();
	    
	    // Add an edge between the two leaves (direction is arbitrarily
	    // chosen to be from leaf 1 to 2)
		// NOTE: from Google Guava "You can add an edge whose incident nodes 
		// have not previously been added to the graph. If they're not already present, 
		// they're silently added to the graph"
	    dend.addEdge(leaf1, leaf2);
	    
	    int id = 1; // (negative) IDs for the internal nodes
	    
	    while (!queue.isEmpty()) {
	      
	    		// Pick an edge at random
	    		EndpointPair<BasicTreeNode> edge = SamplingUtils.uniformFromCollection(rand, dend.getTree().edges());
	    		
	    		// Create a new internal node
			// initially set all internal nodes as left children
			BasicTreeNode internal = BasicTreeNode.createInternal(id++, true);

	    		BasicTreeNode nextLeaf = queue.poll();
	    		
	    		// Remove the selected edge
	    		BasicTreeNode first = edge.nodeU();
	    		BasicTreeNode second = edge.nodeV();
	    		dend.removeEdge(first, second);
	    		
	    		// Add edge from internal node to nextLeaf
	    		// Note: the graph is directed and the directed edges always flow from the internal
	    		// node to the leaves
	    		dend.addEdge(internal, nextLeaf);
	    		
	    		// Add edges from internal node to the incident nodes of chosen edge
	    		// (direction of the edge is chosen arbitrarily at this point)
	    		dend.addEdge(internal, first);
	    		dend.addEdge(internal, second);
	    		
	    		// make second node the right child
	    		second.setLeft(false);

	    }

	    ////////////////////////////////////////////////////////////////////////////////
	    // STAGE 2: Transform the unrooted binary tree with n+1 leaves into a rooted
	    // binary tree with n leaves
	    ////////////////////////////////////////////////////////////////////////////////

	    // Approach: given an unrooted binary tree with n+1 leaves, remove the leaf n+1 
	    // and the edge connecting it to the internal node. Label this internal node as the root.

	    // Get the internal node adjacent to the extra leaf
	    Set<BasicTreeNode> adjNodes = dend.getTree().adjacentNodes(extra);
	    if(adjNodes.size() > 1) {
	    	throw new RuntimeException("The extra leaf has more than one adjacent node.");
	    }

	    BasicTreeNode adjInternal = adjNodes.iterator().next();

	    // Remove the extra leaf 
	    // Note: removeNode(node) also removes all edges incident to node
	    dend.getTree().removeNode(extra);

	    // Reorient the dendrogram, so that edge directions follow from the root
	    // towards the leaves. Also, set the left/right properties.
	    reorientDend(dend, adjInternal);

	    // Set the tree
	    this.setTree(dend.getTree());

	    // Set the internal node that was connected to the extra leaf as the root
	    adjInternal.setLeft(null);
	    this.setRoot(adjInternal);
	    
	    //System.out.println("The randomly sampled dendrogram from the uniform prior: \n" + this.prettyPrint());

	}
	
	/**
	 * Given a dendrogram and its root node, reorient the tree so that all edges are correctly
	 * directed down from the root towards the leaves. Also, set the left/right properties of the 
	 * nodes.
	 * 
	 * @param dend
	 * @param root
	 * 
 	 * Results in a dendrogram with directed edges flowing from root to leaves and with left/right 
	 * properties set arbitrarily.
	 */
	void reorientDend(BasicDendrogram dend, BasicTreeNode root) {

		///////////////////
		// Root node case:
		///////////////////
		
		// Get the nodes adjacent to the root
		Set<BasicTreeNode> adjNodes = dend.getTree().adjacentNodes(root);
		Iterator<BasicTreeNode> it = adjNodes.iterator();
		BasicTreeNode child1 = it.next();
		BasicTreeNode child2 = it.next();
		
		//System.out.println("The nodes adjacent to the root:" + adjNodes);
		//System.out.println("Child 1 is: " + child1);
		//System.out.println("Child 2 is: " + child2);
		
		// Get the edges incident to the root
		Set<EndpointPair<BasicTreeNode>> incEdges = dend.getTree().incidentEdges(root);
		
		//System.out.println("The edges incident to the root:" + incEdges);
		
		// Create a list of nodes in the order needed to remove the edges from dend
		List<BasicTreeNode> toRemove = new ArrayList<BasicTreeNode>();

		for(EndpointPair<BasicTreeNode> edge: incEdges) {
			BasicTreeNode first = edge.nodeU();
			BasicTreeNode second = edge.nodeV();
			toRemove.add(first);
			toRemove.add(second);
		}
		
		//System.out.println("The nodes to be removed: " + toRemove);
		
		// Remove the edges
		for(int i=0; i < toRemove.size() - 1; i++) {
			dend.removeEdge(toRemove.get(i), toRemove.get(i+1));
		}

		// Add the edges back, with proper orientation
		dend.addEdge(root, child1);
		dend.addEdge(root, child2);

		
		// Set the second child as the right
		//child1.setLeft(true);
		child2.setLeft(false);
		
		//System.out.println("The root node:" + root);
		//System.out.println("The left child: " + child1);
		//System.out.println("The right child: " + child2);
		
		///////////////////////
		// Internal node case:
		///////////////////////
		
		// If the children are not leaves, then recurse down both branches
		if(!child1.isLeaf()) { reorientInternalNodes(dend, root, child1);}
		if(!child2.isLeaf()) { reorientInternalNodes(dend, root, child2);}
		
		
	}
	
	/**
	 * Recursive method for reorienting the edges between an internal node and its children nodes.
	 * 
	 * We need as an argument the parent node, so we are able to determine which other adjacent nodes
	 * are the children.
	 * 
	 */
	private static void reorientInternalNodes(BasicDendrogram dend, BasicTreeNode parent, BasicTreeNode current) {
		
		// Get the nodes adjacent to the current node
		Set<BasicTreeNode> adjNodes = dend.getTree().adjacentNodes(current);
		
		// Find all children nodes from this set
		List<BasicTreeNode> childList = new ArrayList<BasicTreeNode>(2);
		
		for(BasicTreeNode node: adjNodes) {
			if(node != parent) { childList.add(node);}
		}
		
		BasicTreeNode child1 = childList.get(0);
		BasicTreeNode child2 = childList.get(1);

		// Get the edges incident to the current node
		Set<EndpointPair<BasicTreeNode>> incEdges = dend.getTree().incidentEdges(current);
		
		// Create a list of nodes in the order needed to remove the edges from dend
		List<BasicTreeNode> toRemove = new ArrayList<BasicTreeNode>();

		for(EndpointPair<BasicTreeNode> edge: incEdges) {
			BasicTreeNode first = edge.nodeU();
			BasicTreeNode second = edge.nodeV();
			toRemove.add(first);
			toRemove.add(second);
		}
		
		// Remove all edges
		// Note: this removes the {parent -> current} edge, which we simply add back
		for(int i=0; i < toRemove.size() - 1; i++) {
			dend.removeEdge(toRemove.get(i), toRemove.get(i+1));
		}

		// Add the edges back, with proper orientation
		dend.addEdge(parent, current);
		dend.addEdge(current, child1);
		dend.addEdge(current, child2);

		// Set the second child as the right
		child2.setLeft(false);
		
		// If a child is not a leaf, then recurse down that branch
		if(!child1.isLeaf()) { reorientInternalNodes(dend, current, child1);}
		if(!child2.isLeaf()) { reorientInternalNodes(dend, current, child2);}
	
	}
	
	/**
	 * Extract and remove a specified number of elements from a list of BasicTreeNodes at random
	 * 
	 */
	private List<BasicTreeNode> getRandomElements(List<BasicTreeNode> list, int numItems, Random rand) {
		
		if(numItems > list.size()) {
			throw new IllegalArgumentException("Cannot extract more elements than available in the given list.");
		}
		// Create a new list for storing the selected elements
		List<BasicTreeNode> newList = new ArrayList<BasicTreeNode>();
		
		for(int i=0; i < numItems; i++) {
			
			// Randomly choose an index between 0 and the size of list
			int randomIdx = rand.nextInt(list.size());
			
			// Add the randomly chosen element to the new list
			newList.add(list.get(randomIdx));
			
			// Remove the chosen element from the original list
			list.remove(randomIdx);
		}
		
		return newList;
	}


}
