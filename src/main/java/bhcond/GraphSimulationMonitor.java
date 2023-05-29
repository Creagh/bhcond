package bhcond;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.eclipse.xtext.xbase.lib.Functions.Function1;

import bayonet.distributions.Random;
import bhcond.GraphAdjMatrix.MutableGraphAdjMatrix;
import blang.core.RealVar;
import blang.inits.experiments.tabwriters.TidilySerializable;
import blang.inits.experiments.tabwriters.TidySerializer.Context;
import briefj.run.Results;
import viz.components.MatrixViz;
import viz.components.TreeViz;
import viz.core.Viz;
import xlinear.Matrix;
import xlinear.MatrixOperations;

public class GraphSimulationMonitor implements TidilySerializable {
  final MutableGraphAdjMatrix simulatedGraph;
  final RealVar alpha, beta;
  final RealVar theta; // mixture parameter
  final RealVar rho; // global prob. parameter from Truncate Geometric compoenent
  final Random random = new Random(1);
  final BasicDendrogram dendogram;
  final ModelVariant variant;
  final int visualizationThinning = 1000;
  static int iter = 0;
  
  public GraphSimulationMonitor(int [] vertices, RealVar alpha, RealVar beta, RealVar theta, RealVar rho, BasicDendrogram dendogram, ModelVariant variant) {
    this.simulatedGraph = new MutableGraphAdjMatrix(vertices);
    this.alpha = alpha;
    this.beta = beta;
    this.theta = theta;
    this.rho = rho;
    this.dendogram = dendogram;
    this.variant = variant;
  }
  
  @Override
  public void serialize(Context context) {
    simulatedGraph.sampleGraph(random, alpha.doubleValue(), beta.doubleValue(), theta.doubleValue(), rho.doubleValue(), dendogram, variant);
    
    // save the simulated graph as a list of edges
    int[][] edges = simulatedGraph.convertToEdgeList();
    try {
		saveSampledEdgeList(edges);
	} catch (IOException e1) {
		System.out.println("Error attempting to save the simulated graph as an edge list.");
		e1.printStackTrace();
	}
    
    //compute the edge density as the proportion of present edges from all possible edges 
    double n = simulatedGraph.getNumVertices(); // store as doubles to avoid integer division
    double e = simulatedGraph.getNumEdges(); 
    //System.out.println(n + " " + e);
    double density = (2 * e) / (n * (n-1)); 
    
    context.recurse(density, "statistic", "density");
    
    // compute the mean vertex degree
    int[] vertices = simulatedGraph.getVertices();
    double degTotal = 0.0; // store as doubles to avoid integer division
    
    for(int v : vertices) {
    	degTotal += simulatedGraph.getNeighbours(v).size();
    }
    double meanDeg = degTotal / vertices.length; 
    
    context.recurse(meanDeg, "statistic", "MeanDegree");
    
    // periodically, print out sampled matrix and tree
    if (iter % visualizationThinning == 0) {
      //saveSampledMatrix();
      //saveSampledTree();
    }
    iter++;
  }
  
  private void saveSampledEdgeList(int[][] arr) throws IOException {
	  // write the 2d array to file
	  StringBuilder builder = new StringBuilder();
	  for(int i = 0; i < arr.length; i++) { // rows i
		  for(int j = 0; j < arr[i].length; j++) { // cols j
			  builder.append(arr[i][j]+""); //append to the output string
			  if(j < arr[i].length - 1) //if this is not the last row element
				  builder.append(","); // add comma
		  }
		  builder.append("\n");//append new line at the end of the row
	  }

	  File edgeListFolder = Results.getFileInResultFolder("sampled-edges");
	  edgeListFolder.mkdir();
	  BufferedWriter writer = new BufferedWriter(new FileWriter(new File(edgeListFolder, "" + iter + ".csv")));
	  writer.write(builder.toString()); //save the string representation of the array
	  writer.close();

  }

private void saveSampledTree() {
    Function1<BasicTreeNode, List<BasicTreeNode>>  children = (BasicTreeNode node) -> {
      return new ArrayList<BasicTreeNode>(dendogram.getChildren(node));
    };
    try {
    	TreeViz<BasicTreeNode> treeViz = new TreeViz<BasicTreeNode>(dendogram.getRoot(), children, Viz.fixHeight(100));
    	File treeFolder = Results.getFileInResultFolder("sampled-trees");
    	treeFolder.mkdir();
    	treeViz.output(new File(treeFolder, "" + iter + ".pdf"));
    } catch (Exception e) {
    	System.out.println("Error occured creating treeViz at iteration " + iter);
    	System.out.println(e); 
    	e.printStackTrace(); 
    }
  }
  
  private void saveSampledMatrix() {
    Matrix copy = MatrixOperations.dense(simulatedGraph.getNumVertices(), simulatedGraph.getNumVertices());
    for (int i = 0; i < simulatedGraph.getNumVertices(); i++)
      for (int j = 0; j < simulatedGraph.getNumVertices(); j++) 
        copy.set(i, j, simulatedGraph.adj[i][j]);
    try {
    	MatrixViz viz = new MatrixViz(copy, MatrixViz.greyScale, Viz.fixHeight(100));
    	File matrixFolder = Results.getFileInResultFolder("sampled-matrices");
    	matrixFolder.mkdir();
    	viz.output(new File(matrixFolder, "" + iter + ".pdf"));
    } catch (Exception e) {
    	System.out.println("Error occured creating MatrixViz at iteration " + iter);
    	System.out.println(e);
    	e.printStackTrace(); 
    }
  }
  
}
