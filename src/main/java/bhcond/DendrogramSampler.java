package bhcond;

import java.util.Iterator;
import java.util.List;
import java.util.Set;

import bayonet.distributions.Random;
import blang.core.LogScaleFactor;
import blang.distributions.Generators;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.SampledVariable;
import blang.mcmc.Sampler;


public class DendrogramSampler implements Sampler {

	  @SampledVariable BasicDendrogram var;

	  @ConnectedFactor List<LogScaleFactor> numericFactors;

	  @Override
	  public void execute(Random rand) {
		  
		  double oldDens = logDensity(); // log density of the dendrogram
		  double newDens; // log density of newly proposed dendrogram
		  double alpha; // the acceptance probability
		  boolean bernoulli; // a Bernoulli RV
		  
		  /**
		   *  Sample a new dendrogram and calculate new log density.
		   *  New dendrogram is proposed by sampling a dendrogram
		   *  via a nearest-neighbour interchange type move.
		   */
		  // get the internal nodes (all nodes excluding leaves)
		  Set<BasicTreeNode> internalNodes = var.getInternalNodes();
		  
		  // choose an internal node uniformly at random
		  int index = rand.nextInt(internalNodes.size());
		  Iterator<BasicTreeNode> iter = internalNodes.iterator();
		  for (int i = 0; i < index; i++) {
		      iter.next();
		  }
		  BasicTreeNode node = iter.next();
		  
		  // if the selected node is the root, exit
		  //if (var.getRoot().getVertexID() == node.getVertexID())
		  if (node.isLeftChild() == null) {
			  return;
		  }
		  
		  // choose choose a child uniformly at random
		  Set<BasicTreeNode> children = var.getChildren(node);
		  
		  index = rand.nextInt(children.size());
		  Iterator<BasicTreeNode> it = children.iterator();
		  for (int i = 0; i < index; i++) {
		      it.next();
		  }
		  BasicTreeNode child = it.next();
		  
		  
		  // perform an interchange move to propose a new dendrogram
		  var.interchange(node, child);
		  newDens = logDensity();
		  
		  /**
		   *  Calculate the acceptance probability
		   */
		  alpha = Math.min(1, Math.exp(newDens - oldDens));
		  
		  /**
		   *  Determine whether to accept the new dendrogram
		   *  Flip a weighted coin: heads = accept, tails = reject
		   */
		  bernoulli = Generators.bernoulli(rand, alpha);
		  
		  if(!bernoulli) {
			  // Reject & undo the interchange move
			  // The same function call as before should revert the changes
			  var.interchange(node, child);
		  }
		  
	  }
	  
	  private double logDensity() {
	    double sum = 0.0;
	    for (LogScaleFactor f : numericFactors)
	      sum += f.logDensity();
	    return sum;
	  }

}
