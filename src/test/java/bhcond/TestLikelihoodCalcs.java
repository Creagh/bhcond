package bhcond;

import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

import bhcond.BasicDendrogram;
import bhcond.BasicTreeNode;
import bhcond.GraphAdjMatrix;

import static bhcond.StaticUtils.*;

public class TestLikelihoodCalcs {
	
	public static void main(String[] args) {
		
		// Build a dendrogram with n leaves
		int numLeaves = 3;
		int[] vertices = IntStream.rangeClosed(1, numLeaves).toArray();
		

		// create empty graph
		GraphAdjMatrix graph = new GraphAdjMatrix(vertices);
		BasicDendrogram d = BasicDendrogram.initDend(graph);
		// note: recall that initDend performs a shuffle of vertices
		
		System.out.println("Number of edges = " + graph.getNumEdges());
		graph.printMatrix(graph.adj);

		System.out.println(d);
		System.out.println(d.prettyPrint());
		
		double alpha = 1.0;
		double beta = 1.0;
		double theta = 0.0;
		double rho = 0.0;
		
		double val = calcLikelihood(d, alpha, beta, theta, rho, graph);
		System.out.println("Log-likelihood = " + val);
		
	}

	public static double calcLikelihood(BasicDendrogram dend, double alpha, double beta, double theta, double rho,
			GraphAdjMatrix graph) {
	      
	      // iterate over all internal nodes
	      // calculate logProbability for each internal node
	      // add to current sum
	      
	      double sum = 0.0;
	      
	      // get all internal tree nodes in dend
	      // Note: there are n-1 internal nodes for n vertices
	      Set<BasicTreeNode> internalNodes = dend.getInternalNodes();
	      
	      for(BasicTreeNode r: internalNodes) {
	      
	        int Lr = dend.countLRleaves(r, true); // number of leaves in left subtree rooted at r
	        int Rr = dend.countLRleaves(r, false); // number of leaves in right subtree rooted at r
	        
	        /**
	        * Calculate E_r, the number of edges in graph whose endpoints have r as their LCA in dend
	        */ 
	        
	        // extract the vertex IDs of the leaf nodes from the left and right subtrees rooted at r
	        List<Integer> left = dend.getLRSubtree(r, true);
	        List<Integer> right = dend.getLRSubtree(r, false);
	        // by default, r must be the LCA between the leaf nodes in the sets 'left' and 'right'
	        
	        // count the number of edges in the graph between the left and right sets of vertices
	        int Er = graph.countEdgesBetween(left, right);
	        
	        /**
	        * Compute the marginalized likelihood p(G | D, alpha, beta, theta, rho)
	        */
	         
	        // First, compute the component due to the Beta - Binomial Mixture (log-scale)
	        // Since it appears in all model variations
	        double binPart = Math.log(1 - theta) + betaBinomialLogDensity(alpha, beta, Lr, Rr, Er);
	         
	        if(theta == 0) { // use the base model, without another mixture component
	          sum += binPart;
	        } else { // use the mixture model
	          if(rho == 0) { // the other mixture component comes from the point mass at zero
	            // Compute the likelihood component from the point mass
	            if(Er == 0) { // when E_r = 0 add both mixture components using the LogSumExp trick
	              sum += bayonet.math.NumericalUtils.logAdd(Math.log(theta), binPart);
	            } else { // when E_r > 0 add the component due to the Beta mixture
	              sum += binPart;	
	            }
	          } else { // the other mixture component comes from the Truncated Geometric
	            // Compute the likelihood component due to the Truncated Geometric (log-scale)
	            double geoPart = Math.log(theta) + truncatedGeometricLogDensity(rho, Lr, Rr, Er);
	            // add both mixture components using the LogSumExp trick
	            sum += bayonet.math.NumericalUtils.logAdd(geoPart, binPart);
	          }
	        }
	      }
	      //System.out.println("Log density = " + sum);
	      return sum;
	    }
	
}