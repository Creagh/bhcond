package bhcond;

import java.io.File;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import blang.inits.ConstructorArg;
import blang.inits.DesignatedConstructor;
import briefj.BriefIO;
import briefj.collections.UnorderedPair;

/**
 * A simple, undirected, unweighted graph, represented as a set of edges.
 * 
 * @author Creagh Briercliffe (creagh.briercliffe@gmail.com)
 *
 */
public class GraphEdgeSet {

	private Set<UnorderedPair<Integer, Integer>> edges; // a set of edges
	private int[] vertices; // an array containing all vertices
	// it will be useful to have these vertices readily available when creating the dendrogram
	
	/**
	 * Read in a graph from a CSV file of edges
	 * 
	 * @param file - CSV file with two columns; each row defines an undirected edge between two vertices,
	 * represented by their integer ID values. E.g. "3, 5" denotes an edge between vertex 3 and vertex 5.
	 * @return GraphEdgeSet
	 */
	@DesignatedConstructor
	public static GraphEdgeSet parse(@ConstructorArg(value = "file") File file) {
		
		Set<UnorderedPair<Integer, Integer>> edges = new LinkedHashSet<>();
		Set<Integer> vertexSet = new LinkedHashSet<Integer>();
		
		// iterate through the lines in the CSV file
		for (List<String> line : BriefIO.readLines(file).splitCSV()) {
			int firstVertex = Integer.parseInt(line.get(0));
			int secondVertex = Integer.parseInt(line.get(0));
			// add an edge
			edges.add(UnorderedPair.of(firstVertex, secondVertex));
			// add vertices to the vertexSet
			vertexSet.add(firstVertex);
			vertexSet.add(secondVertex);
		}
		
		// convert the vertexSet to an array 
		// note: requires java-8
		int[] vertices = vertexSet.stream().mapToInt(i -> i).toArray();
		
		return new GraphEdgeSet(edges, vertices);
	}

	public GraphEdgeSet(Set<UnorderedPair<Integer, Integer>> edges, int[] vertices) {
		super();
		this.edges = edges;
		this.vertices = vertices;
	}

	public Set<UnorderedPair<Integer, Integer>> getEdges() {
		return edges;
	}

	public void setEdges(Set<UnorderedPair<Integer, Integer>> edges) {
		this.edges = edges;
	}

	public int[] getVertices() {
		return vertices;
	}

	public void setVertices(int[] vertices) {
		this.vertices = vertices;
	}
	
	/**
	 * Get the set of edges incident to a given vertex
	 * 
	 * @param vertex
	 * @return incEdges
	 */
	public Set<UnorderedPair<Integer, Integer>> getIncidentEdges(int vertex) {
		// TODO: get subset of edges that contain vertex
		
		return null;
	}
	
}
