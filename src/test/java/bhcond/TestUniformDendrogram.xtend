package bhcond

import java.io.File
import org.junit.Test
import blang.validation.ExactInvarianceTest
import blang.validation.Instance

class TestUniformDendrogram {
  
  @Test
  def void invariance() {
    // create a test case
    val graph = GraphAdjMatrix::parseToGraphAdjMatrix(new File("data/terrorists/terrorist_pairs.csv"))
    val model = (new UniformDendrogram.Builder).setDend(BasicDendrogram::initDend(graph)).build

    // an instance is basically just one Geweke test 
    val instance = new Instance(
      model, 
      // hacky test function: an indicator function on if there is a connection from the root to the node 1
      [if ((dend.tree.successors(dend.root).map[toString].toList).contains("1")) 1.0 else 0.0]
    )
    
    //println(model.dend.tree.successors(model.dend.root));
    //println(model.dend.getChildren(model.dend.root));
   	//println(model.dend);
   	//println(model.dend.getLeaves(model.dend.root));
   	
   	val instance2 = new Instance(
   		model,
        // test function: an indicator that the number of leaves in the left subtree = number of leaves in right subtree
        // i.e. an indicator that the tree is "balanced"
   		[if (dend.countLRleaves(dend.root, true) == dend.countLRleaves(dend.root, false)) 1.0 else 0.0]
   	)

    
    // perform the test
    val test = new ExactInvarianceTest => [
      add(instance)
      // more tests can be added here, multiple comparison correction will be automatically performed
      add(instance2)
    ]
    test.check
  }

}