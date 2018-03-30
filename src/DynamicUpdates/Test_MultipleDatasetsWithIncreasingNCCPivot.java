package projects.DynamicUpdates;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;

public class Test_MultipleDatasetsWithIncreasingNCCPivot {
		
	public static void main(String[] args) throws IOException {	
		
		int numberOfDatasets = 1;
		String[] dataset = new String[numberOfDatasets];
		
		dataset[0] = "elliptical_10_2_bigger";
				
		//select dataset
		int[] datasetsToSelect = {0};
				
		//sampling settings
		int burnIn = 0;
		int numOfSamples = 100;
		int thinning = 1;

		//the concentration parameters for the DP to be tested
		int[] nToTest = {500,550,600,650,700,750,800,850,900,950,1000};
		
		//the fixed alpha when we vary N
		int fixedAlpha = 2;
		
		//the number of measurements (100 is usually a good value)
		int numOfMeasurements = 400;
		
		//set this to 'true' if the similarity to the true solution is of interest.
		//if this is set to 'false' the average similarity to all sample partitions is taken.
		boolean compareMPToTrueSolution = false;
		
		//this adds the final plot code for R
		//currently it creates code for the settings firstAlpha = 1, lastAlpha = 10,
		//algorithms = {0,1,2,3,4,5,6} and measures = {0,1,2,3,4,5,6}.
		//So don't use it when these settings are changed.
		boolean addFinalPlotCode = true;

		//Evaluation of different measures with generated partitions from a Dirichlet Process Gaussian Mixture model
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//select here which similarity measures for the mean partition should be tested:
		//0 = partition distance, 1 = rand index, 2 = adjusted rand index, 3 = Fowlkes-Mallows Index, 4 = Jaccard Index
		//5 = Mirkin Metric, 6 = Wallace Coefficient 1, 7 = Wallace Coefficient 2, 8 = van Dongen Measure, 9 = Normalized Mutual Information 
		//-1 = Pick-A-Cluster + CC-Pivot
		int[] algorithms = {0,9,2,3,4,5,-1}; //{0,1,2,3,4,5,6,7,8};
		
		//select here which measures should be used to calculate the similarity of the mean partition solution to the criterium.
		//The same measures are available as for the algorithms (see description above)
		int[] measures = {11,9,2,3,4,10}; //{0,1,2,3,4,5,6,7,8};
		
		//output dir
		String outputDir = "C:\\Users\\Thomas\\Desktop\\";

		//output file name
		String outputFileName = "Multiple_Datsets_Similarity_Measures_Summary";
		
		int numOfNToTest = nToTest.length;
		int numOfDatasetsToTest = datasetsToSelect.length;
		
		double[][][][] similarityToCriterium = new double[measures.length][numOfNToTest][algorithms.length][numOfMeasurements];
		double[][][] meanSimilarityToCriterium = new double[measures.length][numOfNToTest][algorithms.length];
		double[][][] halfBar95CIForSimilarityToCriterium = new double[measures.length][numOfNToTest][algorithms.length];
		double[][][] elapsedTime = new double[numOfNToTest][algorithms.length][numOfMeasurements];
		double[][] meanElapsedTime = new double[numOfNToTest][algorithms.length];
		double[][] halfBar95CIForElapsedTime = new double[numOfNToTest][algorithms.length];
		double[][][] numOfIterations = new double[numOfNToTest][algorithms.length][numOfMeasurements];
		double[][] meanIterations = new double[numOfNToTest][algorithms.length];
		double[][] halfBar95CIForNumOfIterations = new double[numOfNToTest][algorithms.length];
		int[][] meanPartition = new int[algorithms.length][];
		double[] meanNumOfClusters = new double[numOfNToTest];
		String outputString = "", datafile = "";
		double[][] observations = null, allObservations;
		
		for(int dsIndex = 0; dsIndex < numOfDatasetsToTest;dsIndex++) {
		
			//dataset to test
			datafile = "Datasets\\" + dataset[datasetsToSelect[dsIndex]] + ".csv";
	
			allObservations = CSV.read(datafile);
			
			//define true solution
			int[] trueSolution = new int[allObservations.length];
			for(int i = 0;i<trueSolution.length;i++)
					trueSolution[i] = (int) allObservations[i][allObservations[0].length - 1];
			
			//remove class label
			double[][] newObservations = new double[allObservations.length][allObservations[0].length - 1];
			for(int i = 0;i < allObservations.length;i++)
				for(int j = 0;j < allObservations[0].length - 1;j++)
					newObservations[i][j] = allObservations[i][j];
			allObservations = newObservations;
			
			//First, the parameters for the base distribution ...
			double[] priorMu = 		{0,0};
			double priorKappa = 	priorMu.length;		
			double priorNu = 		priorMu.length;	
			double[][] priorPsi = 	{{1,0},
					  				{0,1}};	
			
			for(int N = 0; N < numOfNToTest; N++) {
				
				//copy the relevant individuals for the current condition N
				observations = new double[nToTest[N]][allObservations[0].length];
				for(int i = 0;i < nToTest[N];i++)
					for(int j = 0;j < allObservations[0].length;j++)
						observations[i][j] = allObservations[i][j];
				
				for(int measurement = 0; measurement < numOfMeasurements; measurement++) {
					System.out.println("start measurement " + (measurement + 1) + " with N " + nToTest[N] + " for dataset " + dataset[datasetsToSelect[dsIndex]]);
				
					DPMixture dpm = new DPMixture(new DependentNormalsBase(priorMu, priorPsi, priorKappa, priorNu, observations), fixedAlpha);
					
					int[][] samples = dpm.simulate(burnIn, numOfSamples, thinning);
					
					//calculate mean number of clusters in all partitions
					for(int i = 0; i < samples.length; i++)
						meanNumOfClusters[N] += getNumOfClusters(samples[i]);
					
					PartitionDistribution PD = new PartitionDistribution(samples);
					PD.prepareArrayToKeepIterations(algorithms.length);
					
					//initialize mean partition candidate
					//PD.setInitialMP(PD.combinedPickAClusterCCPivot());
					PD.setInitialMP(PD.pickACluster());
					
					for(int a = 0;a < algorithms.length;a++) {
						System.out.println("calculating mp by algorithm " + algorithms[a] + " ...");
						
						long startTime = System.nanoTime();
						if(algorithms[a] == -1)
							meanPartition[a] = PD.combinedPickAClusterCCPivot();
						else if(algorithms[a] == 0)
							meanPartition[a] = PD.getMeanPartition1();
						else
							meanPartition[a] = PD.getMeanPartitionDynSim(algorithms[a]);				
						elapsedTime[N][a][measurement] += (System.nanoTime() - startTime) / 1000000;
						
						for(int m = 0; m < measures.length; m++)
							if(compareMPToTrueSolution)
								similarityToCriterium[m][N][a][measurement] += PD.getSimilarity(meanPartition[a], trueSolution, measures[m]);
							else
								similarityToCriterium[m][N][a][measurement] += PD.getAverageSimilarity(meanPartition[a], measures[m]);
					}
					
					//keep number of iterations for every algorithm
					int[] numOfIter = PD.getNumOfIterations();
					for(int a = 0;a < algorithms.length;a++)
						numOfIterations[N][a][measurement] = numOfIter[a];
				}
			}
		}
		
		// until here all datasets should be processed
		// next the results should be calculated ...
			
		for(int N = 0; N < numOfNToTest;N++) {
			
			outputString += "Datasets: ";
			for(int ds = 0; ds < numOfDatasetsToTest;ds++)
				outputString += dataset[datasetsToSelect[ds]] + (ds == numOfDatasetsToTest - 1 ? "" : ",");
			outputString += "\n";			
			outputString	+= "N: " + nToTest[N] + "\n"
							+ "Fixed alpha: " + fixedAlpha + "\n"
							+ "BurnIn: " + burnIn + "\n"
							+ "number of samples: " + numOfSamples + "\n"
							+ "thinning: " + thinning + "\n"
							+ "number of measurements: " + numOfMeasurements + "\n"
							+ "criterium: " + (compareMPToTrueSolution ? "similarity to true solution" : "similarity to sample") + "\n\n";
			
			meanNumOfClusters[N] /= (numOfSamples * numOfMeasurements * numOfDatasetsToTest);
			outputString += "mean number of clusters in all sample partitions for N " + nToTest[N] + ": " + meanNumOfClusters[N] + "\n";
			for(int a = 0;a < algorithms.length;a++) {
				for(int m2 = 0; m2 < numOfMeasurements; m2++) {
					elapsedTime[N][a][m2] /= numOfDatasetsToTest;
					for(int m = 0; m < measures.length; m++) {
						similarityToCriterium[m][N][a][m2] /= numOfDatasetsToTest;
						halfBar95CIForSimilarityToCriterium[m][N][a] = RCodeGenerator.getHalfBarLengthOf95CI(similarityToCriterium[m][N][a]);
						meanSimilarityToCriterium[m][N][a] = RCodeGenerator.getMean(similarityToCriterium[m][N][a]);
					}
				}
			}
			
			for(int a = 0;a < algorithms.length;a++) {
				outputString += "\nmean partition algorithm with similarity measure " + algorithms[a] + "\n";
				for(int m = 0; m < measures.length; m++) {								
					outputString += "average similarity to criterium according to measure " + measures[m] + ": " + meanSimilarityToCriterium[m][N][a];				
					boolean highestValue = true, lowestValue = true;
					for(int a2 = 0;a2 < algorithms.length;a2++)
						if(a2 != a && meanSimilarityToCriterium[m][N][a2] >= meanSimilarityToCriterium[m][N][a])
							highestValue = false;
					for(int a2 = 0;a2 < algorithms.length;a2++)
						if(a2 != a && meanSimilarityToCriterium[m][N][a2] <= meanSimilarityToCriterium[m][N][a])
							lowestValue = false;
					outputString += (highestValue ? " +\n" : (lowestValue ? " -\n" : "\n"));				
				}
				meanElapsedTime[N][a] = RCodeGenerator.getMean(elapsedTime[N][a]);
				halfBar95CIForElapsedTime[N][a] = RCodeGenerator.getHalfBarLengthOf95CI(elapsedTime[N][a]);						
				outputString += "Mean elapsed time (in ms) for algorithm with similarity measure " + algorithms[a] + ": " + meanElapsedTime[N][a] + " (+/- " + halfBar95CIForElapsedTime[N][a] + ")\n";
				meanIterations[N][a] = RCodeGenerator.getMean(numOfIterations[N][a]);
				halfBar95CIForNumOfIterations[N][a] = RCodeGenerator.getHalfBarLengthOf95CI(numOfIterations[N][a]);
				outputString += "Mean number if iterations for algorithm with similarity measure " + algorithms[a] + ": " + meanIterations[N][a] + " (+/- " + halfBar95CIForNumOfIterations[N][a] + ")\n";
			}	
			outputString += "\n";
		}
		
		if(addFinalPlotCode) {
			
			outputString += "# plot code\n";
			outputString += "################################################################################\n";
			
			int[] colors = {1,2,3,4,5,6,8};
			int[] pch = {0,1,2,3,4,5,6};
			String[] measureNames = {"Transfer Distance","Norm. Mutual Information","Adjusted Rand Index","Fowlkes-Mellows Index","Jaccard Index","Mirkin Metric","Wallace Coefficient 1"};
			String[] algorithmNames = {"\"LS: Transfer Distance\"","\"LS: Norm. Mutual Information\"","\"LS: Adjusted Rand Index\"","\"LS: Folwkes-Mellows Index\"","\"LS: Jaccard Index\"","\"LS: Mirkin Metric\"","\"Pick-A-Cluster + CC-Pivot\""};
			
			String plotCode = RCodeGenerator.linePlotsFrom2DimArray(meanElapsedTime, halfBar95CIForElapsedTime, RCodeGenerator.asDouble(nToTest), false, colors, pch, "\"N\"", "\"Elapsed time (in ms)\"", "\"Runtime\"",true,"\"Legend of algorithms\"", algorithmNames);		
			plotCode += RCodeGenerator.linePlotsFrom2DimArray(meanIterations, halfBar95CIForNumOfIterations, RCodeGenerator.asDouble(nToTest), false, colors, pch, "\"N\"", "\"Mean number of iterations\"", "\"Iterations\"",false, "",null);
			
			for(int m = 0; m < measures.length; m++)
				plotCode += RCodeGenerator.linePlotsFrom2DimArray(meanSimilarityToCriterium[m], halfBar95CIForSimilarityToCriterium[m], RCodeGenerator.asDouble(nToTest), false, colors, pch, "\"N\"", "\"Similarity\"", "\"" + measureNames[m] + "\"", false, "",null);
									
			plotCode = RCodeGenerator.plotAsMatrix(3, 3, plotCode);
			
			outputString += plotCode;
		}
			
		writeToFile(outputDir + outputFileName + "_" + System.nanoTime() + ".txt", outputString);
		System.out.println("finished");
	}
	
	//Routine is based on code from user "Windle" published here: https://stackoverflow.com/questions/13707223/how-to-write-an-array-to-a-file-java
	public static void writeToFile(String filename, String data) throws IOException{
		  BufferedWriter outputWriter = null;
		  outputWriter = new BufferedWriter(new FileWriter(filename));
		  outputWriter.write(data);
		  outputWriter.flush();  
		  outputWriter.close();  
	}
	
	public static int getNumOfClusters(int[] partition) {
		boolean[] doesExist  = new boolean[partition.length];
		int numOfClusters = 0;
		for(int i = 0; i < partition.length;i++)
			if(!doesExist[partition[i]]) {
				doesExist[partition[i]] = true;
				numOfClusters++;
			}
		return numOfClusters;
	}
}
