package projects.DynamicUpdates;

/**
 * Class for providing basic CSV reading functions.
 * 
 * @author Thomas J. Glassen
 */

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;

public class CSV {
	
	/*
	 * Routine of martinus:
	 * https://stackoverflow.com/questions/453018/number-of-lines-in-a-file-in-java
	 */
	private static int countLines(String filename) throws IOException {
	    InputStream is = new BufferedInputStream(new FileInputStream(filename));
	    try {
	        byte[] c = new byte[1024];
	        int count = 0;
	        int readChars = 0;
	        boolean empty = true;
	        while ((readChars = is.read(c)) != -1) {
	            empty = false;
	            for (int i = 0; i < readChars; ++i) {
	                if (c[i] == '\n') {
	                    ++count;
	                }
	            }
	        }
	        return (count == 0 && !empty) ? 1 : count;
	    } finally {
	        is.close();
	    }
	}
	
	/*
	 * Routine of Blorgbeard & CraigTP:
	 * https://stackoverflow.com/questions/1102891/how-to-check-if-a-string-is-numeric-in-java
	 */
	public static boolean isNumeric(String str)
	{
	  return str.matches("-?\\d+(\\.\\d+)?");  //match a number with optional '-' and decimal.
	}
	
	/*
	 * This routine uses code from Caffe Latte:
	 * https://stackoverflow.com/questions/18033750/read-one-line-of-a-csv-file-in-java
	 */
	public static double[][] read(String filename, String separator) {
		double[][] data = null;   
		try {
			BufferedReader reader =new BufferedReader(new FileReader(filename));
		
			int numOfRows = countLines(filename);	
		      
			String line = "";
			int i = 0;
			
			boolean firstLineIsChecked = false;
		    while((line=reader.readLine())!=null){
		        String [] rawData =line.trim().split(separator);
		        
		        if(!firstLineIsChecked) {
		        	firstLineIsChecked = true;
		        	if(!isNumeric(rawData[0].trim())) {
		        		data = new double[numOfRows-1][rawData.length];
		        		continue;
		        	} else
		        		data = new double[numOfRows][rawData.length];	
		        }
		        
		        for(int j = 0; j < rawData.length; j++)		        	
		        	if(isNumeric(rawData[j].trim()))
		        			data[i][j] = Double.parseDouble(rawData[j]);
		        	else
		        		data[i][j] = Double.POSITIVE_INFINITY;
		        i++;
		    }
		    
		    reader.close();
        
		}catch(Exception e) {
			System.out.println("Error while trying to read CSV file");
		}
		
        return data;
	}
	
	public static double[][] read(String filename) {
		return read(filename, ",");
	}
}
