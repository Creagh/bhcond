/**
 * 
 */
package bhcond;

/**
 * The nodes of a rooted binary tree (dendrogram)
 * 
 * FOR USE WITH GraphAdjMatrix class.
 * 
 * This object does _not_ store counts like E_r, for each internal node, r.
 * 
 * @author Creagh Briercliffe (creagh.briercliffe@gmail.com)
 *
 */
public class BasicTreeNode {
	
	// Leaf nodes will have an integer identifier;
	// Internal nodes will have a negative integer ID
	private final int id; // vertex ID for a leaf node
	
	private Boolean left; // TRUE if this node is the left child of its parent
	// FALSE if this node is the right child of its parent
	// For the root node, set left to null
	
	/**
	 * Constructor for internal nodes.
	 * 
	 * @param left
	 * @param id - a positive integer, which will be made negative upon instantiation 
	 * @return new BasicTreeNode
	 */
	public static BasicTreeNode createInternal(int id, Boolean left) {
	  return new BasicTreeNode(-id, left);
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + id;
		return result;
	}

	// Compare ID values for checking equality
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		BasicTreeNode other = (BasicTreeNode) obj;
		if (id != other.id)
			return false;
		return true;
	}
	
	@Override
	public String toString() { return "" + id; }

	/**
	 * Constructor for leaf nodes.
	 * 
	 * @param vertexId
	 * @param left
	 * @return new BasicTreeNode
	 */
	public static BasicTreeNode createLeaf(int vertexId, Boolean left) {
    return new BasicTreeNode(vertexId, left);
  }
	
	/**
	 * The internal constructor method called by the specific node constructor (internal or leaf)
	 * 
	 * @param id
	 * @param left
	 */
	private BasicTreeNode(int id, Boolean left) {
		this.id = id;
		this.left = left;
  }


	/**
	 * Get the vertex ID of a node
	 * 
	 */
	public int getVertexID() {
		return(this.id);
	}
	
	/**
	 * Set the left property for the node
	 * 
	 * @param Boolean left - TRUE if left child, FALSE if right, NULL if root.
	 */
	public void setLeft(Boolean left) {
		this.left = left;
	}
	
	/**
	 * Determine if this node is the left child of its parent
	 * 
	 * @return Boolean - TRUE if left, FALSE if right, NULL if root.
	 */
	public Boolean isLeftChild() { return left;}
	
	
	/**
	 * Determine if this node is a leaf;
	 * check that the identifier is not greater than -1.
	 * 
	 * @return boolean
	 */
	public boolean isLeaf() { return id > -1; }
	
}
