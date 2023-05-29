/**
 * 
 */
package bhcond;

import java.util.Set;

import briefj.collections.UnorderedPair;

/**
 * The nodes of a rooted binary tree (dendrogram)
 * 
 * @author Creagh Briercliffe (creagh.briercliffe@gmail.com)
 *
 */
public class TreeNode {
	
	// Let r denote this node:
	private Integer e; // number of edges in G whose endpoints have r as their lowest common ancestor in D (E_r)
	// internal nodes will have values for e; null for leaf nodes.
	private int numLeaves; // total number of leaves in both the left and right subtrees rooted at r (L_r + R_r)
	// leaf nodes will have numLeaves = 1
	
	// For each internal node, r, let leaves(r) denote the set of leaf nodes descended from r.
	// We store the set of edges in G incident to any vertex in leaves(r) as incEdges.
	private Set<UnorderedPair<Integer, Integer>> incEdges;
	
	// Leaf nodes will have an identifier 
	// set to -1 for internal nodes
	private final int id; // vertex ID for a leaf node
	
	/**
	 * Constructor for internal nodes.
	 * 
	 * @param e
	 * @param numLeaves
	 * @param incEdges
	 * @return new TreeNode
	 */
	public static TreeNode createInternal(int e, int numLeaves, Set<UnorderedPair<Integer, Integer>> incEdges)
	{
	  return new TreeNode(e, numLeaves, incEdges, -1);
	}
	
	/**
	 * Constructor for leaf nodes.
	 * 
	 * @param incEdges
	 * @param vertexId
	 * @return new TreeNode
	 */
	public static TreeNode createLeaf(Set<UnorderedPair<Integer, Integer>> incEdges, int vertexId)
  {
    return new TreeNode(null, 1, incEdges, vertexId);
  }
	
	/**
	 * The internal constructor method called by the specific node constructor (internal or leaf)
	 * 
	 * @param e
	 * @param numLeaves
	 * @param incEdges
	 * @param id
	 */
	private TreeNode(Integer e, int numLeaves, Set<UnorderedPair<Integer, Integer>> incEdges, int id) 
	{
    this.e = e;
    this.numLeaves = numLeaves;
    this.incEdges = incEdges;
    this.id = id;
  }

  /**
	 * Get the edge count (E_r) for this node
	 * 
	 * @return integer e
	 */
	public int getEdgeCount() {
		return e;
	}
	
	/**
	 * Get the number of leaves in the subtree rooted from this node
	 * 
	 * @return integer numLeaves
	 */
	public int getNumLeaves() {
		return numLeaves;
	}
	
	
	/**
	 * Get the set of edges incident to any leaf in the subtree rooted at this node
	 * 
	 * @return incEdges
	 */
	public Set<UnorderedPair<Integer, Integer>> getIncEdges() {
		return incEdges;
	}

	/**
	 * Get the vertex ID of the leaf node
	 * 
	 * @return name
	 */
	public int getName() {
		return(this.id);
	}
	
	/**
	 * Set the edge count (E_r) for this node
	 * 
	 * @param value
	 */
	public void setEdgeCount(int value) {
		
		if(isLeaf()) 
			throw new RuntimeException();
	  this.e = value;
	}
	
	/**
	 * Determine if this node is a leaf;
	 * check that the identifier is not -1.
	 * 
	 * @return boolean
	 */
	public boolean isLeaf() { return id != -1; }
	
	/**
	 * Set the number of leaves for the subtree rooted at this node
	 * 
	 * @param value
	 */
	public void setNumLeaves(int value) {
		this.numLeaves = value;
	}

	/**
	 * Set the set of edges incident to any leaf in the subtree rooted at this node
	 * 
	 * @param incEdges
	 */
	public void setIncEdges(Set<UnorderedPair<Integer, Integer>> incEdges) {
		this.incEdges = incEdges;
	}
	
	/** 
	 * Update the incident edge set I_r, for this internal node r, given its children
	 * 
	 * @param TreeNode child1, child2
	 */
	public void updateIncEdges(TreeNode child1, TreeNode child2) {
		
		// Let r = this node, and let s and t be the children of node.
		// I_r = I_s union I_t
		
		// get the sets I_s and I_t of incident edges for the subtrees rooted at t and u, respectively
		Set<UnorderedPair<Integer, Integer>> sIncEdges = child1.getIncEdges();
		Set<UnorderedPair<Integer, Integer>> tIncEdges = child2.getIncEdges();
		
		// compute the union
		sIncEdges.addAll(tIncEdges);
		// update I_r to this union
		setIncEdges(sIncEdges);
	}
	
	/**
	 * Update the count E_r, for this internal node r, given its children
	 * 
	 * Recall that E_r is the number of edges in G whose endpoints
	 * have r as their lowest common ancestor in D.
	 * 
	 * @param TreeNode child1, child2
	 */
	public void updateEdgeCount(TreeNode child1, TreeNode child2) {

		// Let r = this node, and let s and t be the children of node.
		// E_r = | I_s intersect I_t |

		// get the sets I_s and I_t of incident edges for the subtrees rooted at s and t, respectively
		Set<UnorderedPair<Integer, Integer>> sIncEdges = child1.getIncEdges();
		Set<UnorderedPair<Integer, Integer>> tIncEdges = child2.getIncEdges();

		// compute the intersection
		sIncEdges.retainAll(tIncEdges);
		// count the number of elements
		setEdgeCount(sIncEdges.size());
	}
	
}
