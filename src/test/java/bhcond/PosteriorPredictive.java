package bhcond;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Stack;

import au.com.bytecode.opencsv.CSVReader;
import bhcond.BasicDendrogram;


public class PosteriorPredictive {

	public static void main(String[] args) {
		
		// Read in from file: alpha & beta samples
		ArrayList<Double> alpha = parseToArrayList(new File("saved_results/terrorist/___/samples/alpha.csv"));
		ArrayList<Double> beta = parseToArrayList(new File("saved_results/terrorist/___/samples/beta.csv"));
		//System.out.println(alpha);
		System.out.println(beta);
		
		// Read in from file dendrogram samples
		BasicDendrogram dend = parseToDend(new File("saved_results/terrorist/2020-06-09-02-28-01-kJsqCHDE.exec/samples/dend.csv"));
		
		
		// For every triple {alpha, beta, dend}, sample a network from the model
		long seed = 2020L;
		Random rand = new Random(seed);
		
			// For every pair {alpha, beta}, sample an array of internal node probabilities 
		
		// Save the network to file
		// OR calculate statistics here and save stats to file


		
		
	}
	
	
	/**
	 * Read in the sampled dendrograms and store as a BasicDendrogram
	 * 
	 * @param File file - CSV file containing all sampled dendrograms in Newick format, with header
	 * @return BasicDendrogram
	 */
	private static BasicDendrogram parseToDend(File file) {
		BasicDendrogram dend = new BasicDendrogram();
		
		CSVReader reader = null;
		
		try {
			reader = new CSVReader(new FileReader(file));
			String [] nextLine;
			
			// remove header line
			reader.readNext();
			
			// TODO:
			// for now, just grab first dendrogram
			dend = convertNewickToDend(reader.readNext()[1]);
			//while ((nextLine = reader.readNext()) != null) {
				// nextLine[] is an array of values from the line
				// grab only the second column's value
				//list.add(Double.parseDouble(nextLine[1]));
			//}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
		// check:
		//System.out.println(dend.toString());
		
		return dend;
	}
	

	/**
	 * Convert the dendrogram from Newick format to a BasicDendrogram object
	 * 
	 * @param string
	 * @return BasicDendrogram
	 */
	private static BasicDendrogram convertNewickToDend(String string) {
		BasicDendrogram dend = new BasicDendrogram();
		
		// test string
		//string = "(((10,2),317),498765);";
		//string = "((1,2),(3,4));";
		
		Stack<BasicTreeNode> stack = new Stack<BasicTreeNode>(); 
		int id = 1; // ID's for internal nodes
		
		for(int i = 0, n = string.length() ; i < n ; i++) { 
		    char c = string.charAt(i); 
		    
		    if(c == '(') {
		    	// create an internal node & add to stack
		    	BasicTreeNode node = BasicTreeNode.createInternal(id, true);
		    	stack.push(node);
		    	id++;
		    	
		    } else if(Character.isDigit(c)) { // c is a number
		    	
		    	// need to check if c is the first digit in a multi-digit number >9
		    	int j = i;
		    	while(Character.isDigit(string.charAt(j+1))) {
		    		// the vertex ID is > 9
		    		j++;
		    	}
		    	
		    	int vertID = 0;
		    	if(j > i) { // the node ID is a multi-digit number
		    		String sID = string.substring(i, j+1);
		    		vertID = Integer.parseInt(sID);
		    		i = j; // progress the iterator to the end of the number
		    	} else {
		    		vertID = Character.getNumericValue(c);
		    	}
		    	
		    	// create a leaf node & add to stack
		    	BasicTreeNode node = BasicTreeNode.createLeaf(vertID, true);
		    	stack.push(node);
		    	
		    } else if(c == ';') {
		    	// pop the top node from stack
		    	BasicTreeNode node = stack.pop();
		    	// check that the stack is now empty
		    	if(!stack.isEmpty()) { throw new RuntimeException("The stack is not empty."); }
		    	// top node becomes root
		    	node.setLeft(null); // left property null for root node
		    	dend.setRoot(node);
		   
		    } else if(c == ')') {
		    		//connect top 2 nodes from stack to the 3rd from top node by adding edges to dend
		    		BasicTreeNode node1 = stack.pop();
		    		BasicTreeNode node2 = stack.pop();
		    		BasicTreeNode node3 = stack.pop();
		    		
		    		dend.addEdge(node3, node2);
		    		dend.addEdge(node3, node1);
		    		node1.setLeft(false); // make node1 the right child
		    		
		    		// add the internal node back to stack
		    		stack.push(node3);		
		   
		    }
		    // note: if c == ',' then do nothing

		}
		
		return dend;
	}


	/**
	 * Read in the sampled values of alpha/beta and store in an ArrayList
	 * 
	 * @param File file - CSV file containing all sampled values of alpha/beta in the second column, with header
	 * @return Arraylist<Double>
	 */
	private static ArrayList<Double> parseToArrayList(File file) {
		ArrayList<Double> list = new ArrayList<Double>();
		CSVReader reader = null;
		
		try {
			reader = new CSVReader(new FileReader(file));
			String [] nextLine;
			
			// remove header line
			reader.readNext();
			
			while ((nextLine = reader.readNext()) != null) {
				// nextLine[] is an array of values from the line
				// grab only the second column's value
				list.add(Double.parseDouble(nextLine[1]));
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

		return list;
	}
}
