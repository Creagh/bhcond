package bhcond;

import bhcond.BasicDendrogram;
import bhcond.BasicTreeNode;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import org.eclipse.xtext.xbase.lib.Pair;

import com.google.common.collect.Lists;
import com.google.common.graph.EndpointPair;

import bayonet.distributions.Random;
import bayonet.math.SamplingUtils;
import blang.distributions.Generators;

public class TestForwardSampling {
	
	public static void main(String[] args) {
		
		// INPUT PARAMETERS:
		//long seed = 2028L; // seed determining dendrogram and network sampled
		//Random rand = new Random(seed);
		
		//long seed_probs = 2020L; // a seed to fix the internal node probabilities, separate from the dend
		//Random rand_probs = new Random(seed_probs);
		
		int n = 100; // number of vertices (nodes)
		
		//double shape1 = 2.0; // shape & rate of Gamma hyperprior for alpha hyperparameter
		//double rate1 = 10;
		
		//double shape2 = 2.0; // shape & rate of Gamma hyperprior for beta hyperparameter
		//double rate2 = 10;
		
		// HYPERPRIORS:
		//double alpha = Generators.gamma(rand, shape1, rate1);
	    //double beta = Generators.gamma(rand, shape2, rate2);
		
		// Uncomment the above & comment the below 2 lines for random alpha & beta
		
		//double alpha = 0.1;
		//double beta = 1.0;
		
		
		// Parameters for mass run
		int num_sims = 1;
		double[] alpha_vals = new double[] {15.0};
		double[] beta_vals = new double[] {5.0};
		
		long[] seeds = new long[num_sims];

		for(int j=0; j<alpha_vals.length; j++) {
			for(int k=0; k<beta_vals.length; k++) {

				// PRIORS:
				// prior for internal node probabilities
				//double[] p = new double[n-1]; // array containing internal node probabilities
				//for(int r = 0; r < n-1; r++) { p[r] = Generators.beta(rand_probs, alpha_vals[j], beta_vals[k]); }

				System.out.println("alpha = " + alpha_vals[j] + ", beta = " + beta_vals[k]);
				// loop through seeds
				for(int i=0; i<seeds.length; i++) {
					seeds[i] = 2021 + i;
					System.out.println("Seed: " + seeds[i]);
					Random rand = new Random(seeds[i]);
					
					// prior for internal node probabilities
					double[] p = new double[n-1]; // array containing internal node probabilities
					for(int r = 0; r < n-1; r++) { p[r] = Generators.beta(rand, alpha_vals[j], beta_vals[k]); }

					// uniform prior for dendrogram
					BasicDendrogram dend = uniformDend(rand, n);

					// SIMULATE DATA:
					// sample network and store as adjacency matrix
					Pair<int[][], double[][]> pair = sampleNetwork(rand, dend, p);
					int[][] adj = pair.getKey();
					//double[][] pr = pair.getArray2();

					// print the results
					//printValues(seeds[i], n, shape1, rate1, shape2, rate2, alpha, beta, p, dend, adj, pr);

					// write the adjacency matrix to file
					// write the array of internal node probabilities to file
					try {
						saveToFile("data/simulated/toy/", adj, p, n, seeds[i], alpha_vals[j], beta_vals[k]);
					} catch (IOException e) {
						// catch block
						e.printStackTrace();
					}
				}


			}
		}
		
		System.out.println("Complete.");

	    
	}
	
	/**
	 * Write and save the adjacency matrix to a CSV file
	 * Write and save the array of internal node probabilities to another CSV file
	 * 
	 * @param int[][] adj
	 * @param double[] p
	 * @throws IOException 
	 */
	private static void saveToFile(String filename, int[][] adj, double[] p, int n, long seed, double alpha,
			double beta) throws IOException {

		//String date = new SimpleDateFormat("yyyy-MM-dd_hh-mm-ss").format(new Date());
		
		// write the matrix to file
		StringBuilder builder = new StringBuilder();
		for(int i = 0; i < adj.length; i++) { // rows i
		   for(int j = 0; j < adj.length; j++) { // cols j
		      builder.append(adj[i][j]+""); //append to the output string
		      if(j < adj.length - 1) //if this is not the last row element
		         builder.append(","); // add comma
		   }
		   builder.append("\n");//append new line at the end of the row
		}
		
		//BufferedWriter writer = new BufferedWriter(new FileWriter(filename + "n" + n + "seed" + seed + "alpha" +
		//alpha + "beta" + beta + ".csv"));
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename + "n" + n + "seed" + seed + ".csv"));
		writer.write(builder.toString()); //save the string representation of the matrix
		writer.close();
		
		// write the array to file
		builder = new StringBuilder();
		for(int i = 0; i < p.length; i++) { 
		    builder.append(p[i]+"");
			if(i < p.length - 1) { //if this is not the last row element
		         builder.append(",\n");
			}
		}
		builder.append("\n");
		
		//writer = new BufferedWriter(new FileWriter(filename + "probs_n" + n + "seed" + seed + "alpha" +
		//		alpha + "beta" + beta + ".csv"));
		writer = new BufferedWriter(new FileWriter(filename + "probs_n" + n + "seed" + seed + ".csv"));
		writer.write(builder.toString()); //save the string representation of the matrix
		writer.close();
	}

	/** 
	 * Print the values of the input, parameters and simulated data
	 */
	private static void printValues(long seed, int n, double shape1, double rate1, double shape2, double rate2,
			double alpha, double beta, double[] p, BasicDendrogram dend, int[][]adj, double[][] pr) {
		System.out.println("=== Input ===");
		System.out.println("Seed: " + seed);
		System.out.println("Number of nodes: " + n);
		System.out.println("Hyperpriors: \n" + "alpha ~ Gamma(" + shape1 + ", " + rate1 + ")");
		System.out.println("beta ~ Gamma(" + shape2 + ", " + rate2 + ")\n");
		
		System.out.println("=== Hyperparameters ===  \n alpha: " + alpha + "\n beta: " + beta + "\n");
		
		System.out.println("=== Parameters === \n Internal node probabilities: " + Arrays.toString(p) + "\n");
		printDoubleMatrix(pr);
		
		//System.out.print("\n Dendrogram from the uniform prior: \n" + dend.prettyPrint());
		//System.out.println(dend + "\n");
		
		System.out.println("=== Simulated Data ===");
	    System.out.println("Adjacency matrix:");
	    //printMatrix(adj);
	}

	/**
	 * Sample a dendrogram uniformly at random
	 * 
	 * @param rand
	 * @param int n - number of vertices/nodes/leaves
	 * @return
	 */
	private static BasicDendrogram uniformDend(Random rand, int n) {
		// adapted from the code in sampleUniform
		
		// Create a new dendrogram object, consisting of a tree
		BasicDendrogram dend = new BasicDendrogram();

		// Create a list to store the tree nodes
		List<BasicTreeNode> nodes = new ArrayList<BasicTreeNode>();

		// Create a BasicTreeNode for each vertex, and add the nodes to the list
		for(int i=1; i<=n; i++) {
			Boolean left = true; // initially set all leaf nodes as left children
			// create the leaf node
			BasicTreeNode leaf = BasicTreeNode.createLeaf(i, left);
			// add the leaf node to the 	list
			nodes.add(leaf);
		}

		// Add an extra leaf, which will get pruned during the transformation
		// back into a rooted tree
		BasicTreeNode extra = BasicTreeNode.createLeaf(n+1, true);
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
		if (queue.isEmpty()) {return dend;} // no nodes

		BasicTreeNode leaf1 = queue.poll();

		if (queue.isEmpty()) {  // only 1 node
			dend.addNode(leaf1);
			return dend;
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
		dend.reorientDend(dend, adjInternal);

		// Set the tree
		dend.setTree(dend.getTree());

		// Set the internal node that was connected to the extra leaf as the root
		adjInternal.setLeft(null);
		dend.setRoot(adjInternal);

		return dend;
	}

	/**
	 * Sample a network given a dendrogram and a set of internal node probabilties.
	 * Stores the network as a (symmetric) adjacency matrix.
	 * 
	 * @param rand
	 * @param dend
	 * @param double[] p - an array containing the internal node probabilities, in arbitrary order, since they are i.i.d.
	 * @return
	 */
	private static Pair<int[][], double[][]> sampleNetwork(Random rand, BasicDendrogram dend, double[] p) {
		
		int n = p.length + 1;
		int[][] adj = new int[n][n];
		
		// create a matrix to store which internal node gets which probability
		double[][] pr = new double[2][n-1];
		
		// get all internal tree nodes in dend
  		// Note: there are n-1 internal nodes for n vertices
  		Set<BasicTreeNode> internalNodes = dend.getInternalNodes();
  		
  		int idx = 0; // the current index for the internal node prob. p_r from the array
  		
  		for(BasicTreeNode r: internalNodes) {
  			
  			// extract the vertex IDs of the leaf nodes from the left and right subtrees rooted at r
  			List<Integer> left = dend.getLRSubtree(r, true);
  			List<Integer> right = dend.getLRSubtree(r, false);
  			
  			// by default, r must be the LCA between the leaf nodes in the sets 'left' and 'right'
  			
  			// for every pair of leaves from the left and right sets, flip a weighted-coin with prob. p_r,
  			// to determine if they will be joined by an edge
  			for(int u: left) {
  				for(int v: right) {
  	  				Boolean edge = Generators.bernoulli(rand, p[idx]);
  	  				int e = edge ? 1 : 0; // convert edge to 0 (false) or 1 (true)
  	  				
  	  				adj[u-1][v-1] = e;
  	  				adj[v-1][u-1] = e; // adj. matrix is symmetric for undirected networks
  				}
  			}
  			
  			// store the probability and the internal node ID together
  			pr[0][idx] = r.getVertexID();
  			pr[1][idx] = p[idx];
  			
  			idx++; // move to next internal node prob.
  			
  		}
  		
		return new Pair<int[][], double[][]>(adj, pr);
	}
	
	/**
	 * Print an adjacency matrix
	 * 
	 * @param 2d integer matrix
	 */
	public static void printMatrix(int[][] matrix) {
		
		for (int i = 0; i < matrix.length; i++) {
		    for (int j = 0; j < matrix[i].length; j++) {
		        System.out.print(matrix[i][j] + " ");
		    }
		    System.out.println();
		}
	}
	
	/**
	 * Print a matrix of doubles
	 * 
	 * @param 2d double matrix
	 */
	public static void printDoubleMatrix(double[][] matrix) {
	    for (int row = 0; row < matrix.length; row++) {
	    		if(row == 0) {System.out.print("Node ID: ");}
	    		else {System.out.print("p_r:     ");}
	        for (int col = 0; col < matrix[row].length; col++) {
	            System.out.printf("%10f  ", matrix[row][col]);
	        }
	        System.out.println();
	    }
	}

}