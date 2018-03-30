package projects.DynamicUpdates;

/**
 * Class for generating R Code for typical plotting tasks.
 * 
 * @author Thomas J. Glassen
 */

public class RCodeGenerator {
	
	public static String countsToRVector(long count) {
		return "c(1:" + count + ")";
	}
	
	private static String createRVector(String[] arr, int numberLimit) {
		String retStr = "c(";
		int maxNumbersPerLine = 100, numbersCounter = 0;
		for(int i = 0;i < arr.length;i++) {
			numbersCounter++;
			retStr += arr[i];
			if(numberLimit != -1 && i == numberLimit)
				break;
			if(i < arr.length - 1) {
				retStr += ",";
				if(numbersCounter == maxNumbersPerLine) {
					retStr += "\n";
					numbersCounter = 0;
				}
			}
		}
		retStr += ")";
		
		return(retStr);
	}
	
	public static String arrayToRVector(String[] arr) {
		return createRVector(arr,-1);
	}
	
	public static String arrayToRVector(double[] arr) {
		String[] tmpArr = new String[arr.length];
		for(int i = 0; i < arr.length;i++)
			tmpArr[i] = "" + arr[i];
		
		return createRVector(tmpArr,-1);
	}
	
	public static String assignToRVar(String varName, String code) {
		return varName + " = " + code + "\n";
	}
	
	public static String twoDimArrayToRMatrix(double[][] arr) {
		double[] tmp_arr = new double[arr.length * arr[0].length];
		for(int i=0;i< arr.length;i++)
			for(int j=0;j< arr[0].length;j++)
				tmp_arr[i * arr[0].length + j] = arr[i][j];
		
		String retStr = "matrix(" + arrayToRVector(tmp_arr) + ", nrow = " + arr.length + ", byrow = T)";
		
		return retStr;
	}
	
	public static String plotAsMatrix(int rows, int columns, String plotCode) {
		String retStr = "\noldpar = par(mfrow=c(" + rows + "," + columns + "))\n";
		retStr += plotCode;
		retStr += "\npar(oldpar)\n";
		
		return retStr;			
	}
	
	/*
	 * creates R Code for scatter plots for every pair of data line from a 2 dim array, where the datapoints 
	 * for each line are either stored in the rows or the columns of the provided 2 dim array.
	 * 
	 * parameters:
	 * arr:						the 2 dim array, where the datapoints for each line are stored
	 * rowWise:					true if the datapoints of every line is found within every row, if they are found in every column instead, set this to false
	 * colors:					the colors (defined as R color numbers) for each line
	 * xlab:					the label of the x-axis
	 * ylab:					the label of the y-axis
	 * main:					the title of the plot
	 */
	public static String scatterPlotsFrom2DimArray(double[][] arr, boolean rowWise, int[] colors, String xlab, String ylab, String main) {
		
		String retStr = "\n#Scatter plots\n";
		retStr += "####################################################################\n\n";
		
		retStr += assignToRVar("data", twoDimArrayToRMatrix(arr));
		retStr += assignToRVar("colors", arrayToRVector(asDouble(colors)));
		
		retStr += "for(i in 1:" + arr.length + "){\n"; 
		retStr += "    for(j in i:" + arr[0].length + "){\n";
		retStr += "        plot(data[" + (rowWise ? "i" : "") + "," + (rowWise ? "" : "i") + "], data[" + (rowWise ? "j" : "") + "," + (rowWise ? "" : "j") + "], xlab = " 
				+ xlab + ", ylab = " + ylab + ", main = " + main + ", col = colors)\n";
		retStr += "    }\n";
		retStr += "}\n\n";
		
		retStr += "####################################################################\n";
		
		return retStr;
	}
	
	/*
	 * creates R Code for a line plot (multiple lines in one plot), where the datapoints 
	 * for each line are either stored in the rows or the columns of a provided 2 dim array.
	 * 
	 * parameters:
	 * data:					the 2 dim array, where the datapoints for each line are stored
	 * dev:						a 2 dim array, where the +/- deviations (e.g. half bar of 95% CI) of each datapoint is stored
	 * xValues:					the values on the x-axis for each datapoint
	 * rowWise:					true if the datapoints of every line is found within every row, if they are found in every column instead, set this to false
	 * colors:					the colors (defined as R color numbers) for each line
	 * pch:						the symbol (defined as R symbol numbers) for each point
	 * xlab:					the label of the x-axis
	 * ylab:					the label of the y-axis
	 * main:					the title of the plot
	 * addPrecedingLegendPlot:	plots first an empty plot that contains only a legend if this is set to true, nothing is precedingly plotted if this is set to false
	 * legendTitle:				title of th legend if addPrecedingLegendPlot = true, if addPrecedingLegendPlot = false this parameter is ignored (assign "" then)
	 * legendLineNames:			names of the lines in the legend if addPrecedingLegendPlot = true, if addPrecedingLegendPlot = false this parameter is ignored (assign 'null' then)
	 */
	public static String linePlotsFrom2DimArray(double[][] data, double[][] dev, double[] xValues, boolean rowWise, int[] colors, int[] pch, String xlab, String ylab, String main, boolean addPrecedingLegendPlot, String legendTitle, String[] legendLineNames) {
		
		String retStr = "\n#Line plot\n";
		retStr += "####################################################################\n\n";
		
		retStr += assignToRVar("data", twoDimArrayToRMatrix(data));
		retStr += assignToRVar("dev", twoDimArrayToRMatrix(dev));
		retStr += assignToRVar("colors", arrayToRVector(asDouble(colors)));
		retStr += assignToRVar("pchs", arrayToRVector(asDouble(pch)));
		retStr += assignToRVar("xValues", arrayToRVector(xValues));
		
		retStr += "yValues = c(";
		int numOfLines = (rowWise ? data.length : data[0].length);
		for(int i = 0; i < numOfLines;i++)
			retStr += "data[" + (rowWise ? i + 1 : "") + "," + (rowWise ? "" : i + 1) + "]" + (i < numOfLines - 1 ? "," : "");
		retStr += ")\n\n";
		
		if(addPrecedingLegendPlot) {
			retStr += "plot(rep(xValues," + numOfLines + "), yValues, " +
					 "type=\"n\", axes=FALSE, xlab=\"\", ylab=\"\", main = " + legendTitle + ", col = 0)\n";
			retStr += "legend(\"center\", legend=" + arrayToRVector(legendLineNames) + ", col=c(1:6,8), pch=0:6, lty = 1, cex = 0.6, pt.cex = 1)\n\n";
		}
		
		retStr += "plot(rep(xValues," + numOfLines + "), yValues, xlab = " 
				+ xlab + ", ylab = " + ylab + ", main = " + main + ", col = 0)\n";
		
		retStr += "for(line in 1:" + numOfLines + "){\n";
		retStr += "    points(xValues, data[" + (rowWise ? "line" : "") + "," + (rowWise ? "" : "line") + "], pch = pchs[line], col = colors[line], cex = 0.8)\n";
		retStr += "    lines(xValues, data[" + (rowWise ? "line" : "") + "," + (rowWise ? "" : "line") + "], col = colors[line])\n";
		retStr += "    #Error bars technique by user 'Laryx Decidua': https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars\n";
		retStr += "    arrows(xValues, data[" + (rowWise ? "line" : "") + "," + (rowWise ? "" : "line") + "]-dev[" + (rowWise ? "line" : "") + "," + (rowWise ? "" : "line") + "], xValues, " +
					"data[" + (rowWise ? "line" : "") + "," + (rowWise ? "" : "line") + "]+dev[" + (rowWise ? "line" : "") + "," + (rowWise ? "" : "line") + "], length=0.025, angle=90, code=3, col = colors[line])\n";
		retStr += "}\n\n";
		
		retStr += "####################################################################\n";
		
		return retStr;
	}
	
	/*
	 * This version calculates the error bars by itself.
	 * It's the same as the routine with parameters data and dev, but this time with all measurements.
	 * The last dimension of 'measurements' keeps the repeated measurements for the conditions in the first
	 * two dimensions.
	 */
	public static String linePlotsFrom2DimArray(double[][][] measurements, double[] xValues, boolean rowWise, int[] colors, int[] pch, String xlab, String ylab, String main, boolean addPrecedingLegendPlot, String legendTitle, String[] legendLineNames) {
		double[][] data = new double[measurements.length][measurements[0].length];
		double[][] dev = new double[measurements.length][measurements[0].length];
		
		for(int i=0;i < measurements.length;i++)
			for(int j=0;j < measurements.length;j++) {
				data[i][j] = getMean(measurements[i][j]);
				dev[i][j] = getHalfBarLengthOf95CI(measurements[i][j]);
			}
		
		return linePlotsFrom2DimArray(data, dev, xValues, rowWise, colors, pch, xlab, ylab, main, addPrecedingLegendPlot, legendTitle, legendLineNames);
	}

	/*
	 * This is a shortcut version which sets standard values for most parameters.
	 */
	public static String linePlotsFrom2DimArray(double[][][] measurements, double[] xValues, boolean rowWise) {
		int numOfLines = (rowWise ? measurements.length : measurements[0].length);
		int[] colors = new int[numOfLines], pch = colors.clone();
		for(int i=0;i<colors.length;i++) {
			colors[i] = i + 1;
			pch[i] = i;
		}
			
		return linePlotsFrom2DimArray(measurements, xValues, rowWise, colors, pch, "\"\"", "\"\"", "\"\"", false, "\"\"", null);
	}
	
	//////////////////////////////////////////////////////////////////////////
	
	public static double getMean(double[] values) {
		double sumOf = 0;
		for(int i = 0;i < values.length;i++)
			sumOf += values[i];
		return sumOf / values.length;
	}
	
	public static double getBesselCorrectedVar(double[] values) {
		double mean = getMean(values);
		double sumOfSquaredDevs = 0;
		for(int i = 0;i < values.length;i++)
			sumOfSquaredDevs += (values[i] - mean) * (values[i] - mean) ;
		return sumOfSquaredDevs / (values.length - 1);
	}
	
	public static double getBesselCorrectedStd(double[] values) {
		return Math.sqrt(getBesselCorrectedVar(values));
	}
	
	//Calculates the length of a half bar of the 95% confidence interval according to:
	//https://de.wikipedia.org/wiki/Konfidenzintervall#.C3.9Cbersicht_f.C3.BCr_stetige_Verteilungen
	public static double getHalfBarLengthOf95CI(double[] values) {
		double ZValueFor5PercentAlpha = 1.959964;
		
		//calculate standard error
		double stdErr = getBesselCorrectedStd(values) / Math.sqrt(values.length);
		
		//calculate half bar length
		double halfBarLength = ZValueFor5PercentAlpha * stdErr;
		
		return halfBarLength;		
	}
	
	//////////////////////////////////////////////////////////////////////////
	
	public static double[] asDouble(int[] arr) {
		double[] tmp_arr = new double[arr.length];
		for(int i=0;i< arr.length;i++)
			tmp_arr[i] = arr[i];
		
		return tmp_arr;
	}
}
