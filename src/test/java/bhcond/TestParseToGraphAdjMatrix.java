package bhcond;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class TestParseToGraphAdjMatrix {
	
	public static void main(String[] args) {
		
		String fileName = "data/example/example.csv";
		File file= new File(fileName);
		
		GraphAdjMatrix graph = GraphAdjMatrix.parseToGraphAdjMatrix(file);
		
		System.out.println("Number of vertices: " + graph.getNumVertices());
		System.out.println("Number of edges: " + graph.getNumEdges());
		System.out.println("Max vertex number: " + graph.getMaxVertexNum());
		
		int[] vertices = graph.getVertices();
		
		for(int i : vertices) {
			System.out.println("vertex " + i);
		}
		
		BasicDendrogram dend = BasicDendrogram.initDend(graph);
		System.out.println(dend.toString());
		
		List<Integer> leaves = dend.getLeaves(dend.getRoot());
		
		for(int i : leaves) {
			System.out.println("leaf ID " + i);
		}
		
		int[][] edges = graph.convertToEdgeList();
		System.out.println(Arrays.deepToString(edges).replace("], ", "]\n"));
	}

}
