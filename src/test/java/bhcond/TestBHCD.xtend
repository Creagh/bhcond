package bhcond

import org.junit.Test
import blang.validation.ExactInvarianceTest
import blang.validation.Instance
import static blang.types.StaticUtils.*
import bhcond.GraphAdjMatrix.MutableGraphAdjMatrix

class TestBHCOND {
	
  @Test
  def void invariance() {
    // create a test case
    val graph = new MutableGraphAdjMatrix(20)
    val alpha = latentReal
    val beta = latentReal
    val theta = latentReal
    val rho = latentReal
    val variant = ModelVariant.BINOM
    val model = (new BHCOND.Builder).setGraph(graph).setDend(BasicDendrogram::initDend(graph)).setAlpha(alpha).setBeta(beta).setTheta(theta).setRho(rho).setVariant(variant).build

    // an instance is basically just one Geweke test 
    val instance = new Instance(
      model, 
      // quick test function: an indicator function on if the number of edges is greater than the number of vertices
      [if (graph.getNumEdges > graph.getNumVertices) 1.0 else 0.0], 
      [it.alpha.doubleValue], 
      [it.beta.doubleValue],
      [(graph.getNumEdges / ((graph.getNumVertices)*(graph.getNumVertices - 1) / 2)).doubleValue]
   	)
   	
    // perform the test
    val test = new ExactInvarianceTest => [
      add(instance)
    ]
    test.nPosteriorSamplesPerIndep = 20
    test.nIndependentSamples = 20_000
    test.check
  }

}