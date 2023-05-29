package bhcond;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

/**
 * File Input taken from https://stackoverflow.com/questions/40074840/reading-a-csv-file-into-a-array
 * @author Michael Lihs
 * 
 * Test Code: read in an input file from CSV and find the number of unique values.
 * Convert the values to a reduced form.
 * 
 * @author Creagh Briercliffe
 *
 */

public class TestFileInput {

	public static void main(String[] args) {
        
		String fileName= "data/grass_web/grass_web_pairs.csv";
        File file= new File(fileName);

        // create a 2-dimensional array of strings
        List<List<String>> lines = new ArrayList<>();
        Scanner inputStream;

        try{
            inputStream = new Scanner(file);

            while(inputStream.hasNext()){
                String line= inputStream.next();
                String[] values = line.split(",");
                // add the currently parsed line to the 2-dimensional string array
                lines.add(Arrays.asList(values));
            }

            inputStream.close();
        }catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // iterate through the 2-dimensional array
        // transfer all values into a single array
        int lineNo = 1;
        List<Integer> values = new ArrayList<Integer>();
        // also, store the columns separately 
        List<Integer> col1 = new ArrayList<Integer>();
        List<Integer> col2 = new ArrayList<Integer>();
        
        for(List<String> line: lines) {
            int columnNo = 1;
            for (String value: line) {
                //System.out.println("Line " + lineNo + " Column " + columnNo + ": " + value);
            		values.add(Integer.parseInt(value));
            		
            		if(columnNo == 1) {
            			col1.add(Integer.parseInt(value));
            		} else {
            			col2.add(Integer.parseInt(value));
            		}
            		
                columnNo++;
            }
            lineNo++;
        }
        
        //System.out.println(values);
        System.out.println(col1);
        System.out.println(col2);
        
        // remove duplicates
        Set<Integer> uniqueValues = new HashSet<Integer>(values);
        	System.out.println("Unique value count: " + uniqueValues.size());
        //System.out.println(uniqueValues);
        	
        	
        	// convert values to reduced form
        	int[] arr = uniqueValues.stream().mapToInt(i->i).toArray();
        	int[] oldArr = Arrays.copyOf(arr, arr.length); // create deep copy of array before converting
        	convert(arr, arr.length);
    	
        	System.out.println(Arrays.toString(oldArr));
    		System.out.println(Arrays.toString(arr));
        
    }
	
    /**
     * Convert the array of values to a reduced form.
     * 
     * Given an array with n distinct elements, convert the given array to a form where
     * all elements are in range from 1 to n. The order of elements is same, i.e., 1
     * is placed in place of smallest element, 2 is placed for second smallest element,
     * … n is placed for largest element.
     * 
     * Source: https://www.geeksforgeeks.org/convert-an-array-to-reduced-form-set-1-simple-and-hashing/
     */
	public static void convert(int[]arr, int n) {
        
        // Create a temp array and copy contents 
        // of arr[] to temp 
        int temp[] = arr.clone(); 
      
        // Sort temp array 
        Arrays.sort(temp); 
      
        // Create a hash table. 
        HashMap<Integer, Integer> umap = new HashMap<>(); 
      
        // One by one, insert elements of sorted 
        // temp[] and assign them values from 1 to n 
        int val = 1; 
        for (int i = 0; i < n; i++) 
        		umap.put(temp[i], val++); 
      
        // Convert array by taking positions from umap 
        for (int i = 0; i < n; i++) {
        		arr[i] = umap.get(arr[i]); 
        }
	}
}
