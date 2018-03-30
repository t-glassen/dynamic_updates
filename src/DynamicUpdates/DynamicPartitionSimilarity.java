package projects.DynamicUpdates;

/**
 * Calculates the partition similarity between two partitions 
 * on various measures and allows for dynamic similarity updates 
 * if an individual is moved from one cluster into another.
 * 
 * @author Thomas J. Glassen
 */

/*
 * This version has a worst case runtime of the unique cluster nr search in O(C), 
 * where C is the highest reached number of clusters in the whole optimization procedure  
 */

public class DynamicPartitionSimilarity {
	
	private int[][] partitions;
	
	private int[] numOfClusters;
	private int[][] clusterSizes;
	private int firstEmptyCluster;
	private int maxExpectedNumOfClusters;
	
	private long[][] ISM; //current Intersection Size Matrix
	private long[] MM; //current Mismatch Matrix
	
	private double[] entropies;
	private double jointEntropy;
	
	public DynamicPartitionSimilarity(int[] partition1,int[] partition2){
		partitions = new int[2][];
		partitions[0] = partition1.clone();
		partitions[1] = partition2.clone();
		
		calcIntersectionSizeMatrix();
		calcMismatchMatrix();
		calcEntropies();
		
		firstEmptyCluster = numOfClusters[1] + 1;
	}
	
	private void increaseSizeOfArrays() {		  
		//set the maximum expected number of clusters. Keep this low to avoid heap size errors
		maxExpectedNumOfClusters = Math.min(Math.max(numOfClusters[0], numOfClusters[1] + 1) * 2, partitions[1].length + 1);
		 
		int[][] newClusterSizes = new int[2][maxExpectedNumOfClusters];
		long[][] newISM = new long[numOfClusters[0]][maxExpectedNumOfClusters];
		  	
	  	for(int j = 0; j < clusterSizes[0].length;j++) {
	  		newClusterSizes[0][j] = clusterSizes[0][j];
	  		newClusterSizes[1][j] = clusterSizes[1][j];
	  		for(int i = 0; i < numOfClusters[0];i++)
	  			newISM[i][j] = ISM[i][j];
	  	}
	  	clusterSizes = newClusterSizes; ISM = newISM;
	}
	
	/**
	 * calculates the unique ascending cluster indices of partition 2
	 * 
	 * @return the ascending cluster indices
	 */
	public int[] getAscendingClusterNrs(){
		int counter = numOfClusters[1];
		int[] clusterNrs = new int[counter];
		for(int i=0;counter > 0;i++)
			if(clusterSizes[1][i] > 0){
				clusterNrs[clusterNrs.length - counter] = i+1;
				counter--;
			}
		return clusterNrs;
	}
	
	/**
	 * Provides the cluster index of an object in partition 2
	 * 
	 * @return the cluster index of obj 
	 */
	public int getClusterNrOfObj(int index){
		return partitions[1][index];
	}
	  
	public int getSizeOfClusterOfObj(int index){
		return clusterSizes[1][partitions[1][index] - 1];
	}
	  
	public int[] getPartition(){
		return partitions[1];
	}
	
	//runtime: O(1)
	public void movObjTo(int index, int k){
		
		//get relevant indices	
		int i = partitions[0][index] - 1;
		int j = partitions[1][index] - 1;
		  
		if(j == k - 1)
			return; //nothing to do		
		
		//update Mismatch Matrix
		MM[0] = MM[0] - (ISM[i][j] - 1) + ISM[i][k-1];
    	MM[1] = MM[1] - (clusterSizes[1][k-1] - ISM[i][k-1]) + (clusterSizes[1][j] - ISM[i][j]);
    	MM[2] = MM[2] - ISM[i][k-1] + (ISM[i][j] - 1);
    	MM[3] = MM[3] - (clusterSizes[1][j] - ISM[i][j]) + (clusterSizes[1][k-1] - ISM[i][k-1]);

    	//update entropies
    	updateEntropy(1,j,k);
    	updateJointEntropy(i,j,k);
    	
    	//update cluster index
    	partitions[1][index] = k;
    	
    	//update cluster sizes
		clusterSizes[1][j]--;
		clusterSizes[1][k-1]++;
		
		//update number of clusters
		if(clusterSizes[1][j] == 0) {
			numOfClusters[1]--;
			//keep cluster nr if j + 1 is lower than firstEmptyCluster
			if(firstEmptyCluster > j + 1)
				firstEmptyCluster = j + 1;
		}
		if(clusterSizes[1][k-1] == 1)
			numOfClusters[1]++;
    	
    	//update Intersection Size Matrix
    	ISM[i][j]--;
    	ISM[i][k-1]++;
	}
	
	//runtime: O(C) (amortized O(1) if it is called once after C movObjTo)
	public void movObjToNew(int index){			  
		//get the next empty cluster. Runtime: O(C)
		int emptyClusterNr = firstEmptyCluster;
		firstEmptyCluster++;
		  
		while(firstEmptyCluster <= clusterSizes[1].length && clusterSizes[1][firstEmptyCluster - 1] > 0)
			firstEmptyCluster++;
		  
		//increase size of arrays if maxExpectedNumOfClusters is nearly reached
		if(firstEmptyCluster == maxExpectedNumOfClusters)
			increaseSizeOfArrays();
		  
		movObjTo(index, emptyClusterNr);
	}
	
	//returns the similarity between two partitions calculated via a specified measure
	public double getSimilarity(int measure) {
		if(measure == 1)		
			return getRandIndex();
		else if(measure == 2)
			return getAdjustedRandIndex();
		else if(measure == 3)
			return getFowlkesMallowsIndex();
		else if(measure == 4)
			return getJaccardIndex();
		else if(measure == 5)
			return getMirkinMetric();
		else if(measure == 6)
			return getWallaceCoefficient1();
		else if(measure == 7)
			return getWallaceCoefficient2();
		else if(measure == 8)
			return getVanDongenMeasure();
		else if(measure == 9)
			return getNormalizedMutualInformation1();
		else if(measure == 10)
			return getBoundedMirkinSimilarity();
		else if(measure == 12)
			return getVariationOfInformation();
		else if(measure == 13)
			return getChiSquaredCoefficient();
		else if(measure == 14)
			return getFMeasure();
		else if(measure == 15)
			return getMeilaHeckermanMeasure();
		else if(measure == 16)
			return getNormalizedMutualInformation2();
		else
			return Double.NaN;
	}

	public String getNameOfMeasure(int measure) {
		String[] measureNames =	{
								"NaN",
								"Rand Index",
								"Adjusted Rand Index",
								"Fowlkes Mallows Index",
								"Jaccard Index",
								"Mirkin Metric",
								"Wallace Coefficient 1",
								"Wallace Coefficient 2",
								"Van Dongen Measure",
								"Normalized Mutual Information 1 (Strehl & Gosh)",
								"Bounded Mirkin Similarity",
								"NaN",
								"Variation of Information",
								"Chi Squared Coefficient",
								"F Measure (Larsen & Aone)",
								"Meila-Heckerman Measure",
								"Normalized Mutual Information 2 (Fred & Jain)"
								};
		return measureNames[measure];
	}
	
	//calculates the Rand Index for two partitions of the same size
	//according to https://en.wikipedia.org/wiki/Rand_index
	//runtime: O(1)
	private double getRandIndex() {	
		return (MM[0] + MM[1]) / (double) (MM[0] + MM[1] + MM[2] + MM[3]);
	}
    
    //calculates the Adjusted Rand Index for two partitions of the same size according to:
    //Santos & Embrechts (On the Use of the Adjusted Rand Index as a Metric for Evaluating Supervised Classification; 2009)
	//runtime: O(1)
	private double getAdjustedRandIndex() {
		long exp = (MM[0] + MM[2]) * (MM[0] + MM[3]) + (MM[3] + MM[1]) * (MM[2] + MM[1]);
		long bin = partitions[0].length * (partitions[0].length - 1) / 2;
		//if we have a one cluster solution in both partitions the denominator becomes zero (but also the nominator)
		//we interpret this as an ARI = 0
		if(bin * bin - exp != 0) {
			return (bin * (MM[0] + MM[1]) - exp) / (double) (bin * bin - exp);
		}else
			return 0;				
    }
	
	//calculates the Fowlkes-Mallows Index for two partitions of the same size
	//according to Wagner & Wagner (Comparing Clusterings; 2007)
    //runtime: O(1)
	private double getFowlkesMallowsIndex() {
		return MM[0] / (Math.sqrt( (MM[0] + MM[2]) * (MM[0] + MM[3]) ) );
	}
	
	//calculates the Jaccard index for two partitions of the same size
	//according to Wagner & Wagner (Comparing Clusterings; 2007)
	//runtime: O(1)
	private double getJaccardIndex() {
		return MM[0] / (double) (MM[0] + MM[2] + MM[3]);
	}
	
    //calculates the Mirkin Metric for two partitions of the same size
  	//according to Wagner & Wagner (Comparing Clusterings; 2007)
	//runtime: O(1)
	private double getMirkinMetric() {
		return 2 * (MM[2] + MM[3]);
    }
	
	//calculates the Bounded Mirkin Metric for two partitions of the same size
	//according to Meila (Comparing Clusterings - An Axiomatic View; 2005)
	//subtracted from 1
	//runtime: O(1)
	private double getBoundedMirkinSimilarity() {
		return 1 - getMirkinMetric() / (partitions[0].length * partitions[0].length);
	}
    
    //calculates the Wallace's coefficient (A -> B) according to:
    //http://www.comparingpartitions.info/index.php?link=Tut8
    //Wallace (A Method for Comparing Two Hierarchical Clusterings: Comment; 1983)
    //runtime: O(1)
	private double getWallaceCoefficient1() {    	
		return MM[0] / (double) (MM[0] + MM[2]);
    }
    
    //calculates the Wallace's coefficient (B -> A) according to:
    //http://www.comparingpartitions.info/index.php?link=Tut8
    //Wallace (A Method for Comparing Two Hierarchical Clusterings: Comment; 1983)
    //runtime: O(1)
	private double getWallaceCoefficient2() {    	
		return MM[0] / (double) (MM[0] + MM[3]);
    }
	
	//calculates the Van Dongen Measure for two partitions of the same size
  	//according to Wagner & Wagner (Comparing Clusterings; 2007)
  	//runtime: O(C^2)
    private double getVanDongenMeasure() {
    	long sumOfMaxI = 0, sumOfMaxJ = 0;
    	long[] maxJ = new long[ISM[0].length];
		for(int i=0;i< ISM.length;i++) {
			long maxI = 0;
			for(int j=0;j< ISM[0].length;j++) {
				if(ISM[i][j] > maxI)
					maxI = ISM[i][j];
				if(ISM[i][j] > maxJ[j])
					maxJ[j] = ISM[i][j];
			}
			sumOfMaxI += maxI;
		}
		for(int j=0;j<maxJ.length;j++)
			sumOfMaxJ += maxJ[j];
		
		return 2 * partitions[0].length - sumOfMaxI - sumOfMaxJ;
    }
    
	//calculates the Chi Squared Coefficient for two partitions of the same size
  	//according to Wagner & Wagner (Comparing Clusterings; 2007)
  	//runtime: O(C^2)
    private double getChiSquaredCoefficient() {
    	double chiSquare = 0;
    	for(int i=0;i< ISM.length;i++) {
    		for(int j=0;j< ISM[0].length;j++) {
    			if(clusterSizes[0][i] != 0 && clusterSizes[1][j] != 0) {
    				double Eij = clusterSizes[0][i] * clusterSizes[1][j] / (double) partitions[0].length;
    				chiSquare += (ISM[i][j] - Eij) * (ISM[i][j] - Eij) / Eij;
    			}
    		}
    	}
    	return chiSquare;
    }
    
    //There is a dissensus between the correct calculation between:
    // Meila (Comparing clusterings—an information based distance; 2007)
    // Wagner & Wagner (Comparing Clusterings; 2007)
    // Zhao & Karypis (Hierarchical Clustering Algorithms for Document Datasets; 2004)
    //The implementation of the two points of dissensus 
    //((1) divide by C or by N at the end, (2) ISM[i][j] or |C_i| * |C'_j| in the nominator of F)
    //in this version is agreed by two of the three mentioned author groups
    //and seems most plausible by reviewing the original article of
    //Larsen & Aone (Fast and Effective Text Mining Using Linear-time Document Clustering, 1999)
    //The variant for (2) of Wagner & Wagner (2007) seems to be wrong, when looking at their own variable definitions.
    //runtime: O(C^2)
    private double getFMeasure() {
    	double F_overall = 0;
    	for(int i=0;i< ISM.length;i++) {
    		double rowFMax = 0;
    		for(int j=0;j< ISM[0].length;j++) { 
    			//Variant of Meila (2007) and Zhao & Karypis (2004)
    			double F = 2 * ISM[i][j] / (double) (clusterSizes[0][i] + clusterSizes[1][j]);
    			//Variant of Wagner & Wagner (2007)
    			//double temp = 2 * clusterSizes[0][i] * clusterSizes[1][j] / (double) (clusterSizes[0][i] + clusterSizes[1][j]);
    			rowFMax = (F > rowFMax ? F : rowFMax);    			
    		}
    		F_overall += clusterSizes[0][i] * rowFMax;
    	}
    	//Variant of Zhao & Karypis (2004) and Wagner & Wagner (2007)
    	return F_overall / (double) partitions[0].length;
    	//Variant of Meila (2007)
    	//return F_overall / (double) numOfClusters[0];
    }
    
	//calculates the Meila-Heckerman Measure for two partitions of the same size
  	//according to Wagner & Wagner (Comparing Clusterings; 2007)
  	//runtime: O(C^2)
    private double getMeilaHeckermanMeasure() {
    	double MHM = 0;
    	for(int i=0;i< ISM.length;i++) {
    		double rowMHMMax = 0;
    		for(int j=0;j< ISM[0].length;j++)
    			rowMHMMax = (ISM[i][j] > rowMHMMax ? ISM[i][j] : rowMHMMax);
    		MHM += rowMHMMax;
    	}
    	return MHM / (double) partitions[0].length;
    }
    
    //calculates the Mutual Information for two partitions of the same size
  	//according to Gosh & Acharya (Cluster ensembles; 2011)
  	//runtime: O(1)
    private double getMutualInformation() {
    	return entropies[0] + entropies[1] - jointEntropy;
    }
    
    //calculates the Strehl & Gosh Normalized Mutual Information for two partitions of the same size
  	//according to Gosh & Acharya (Cluster ensembles; 2011)
  	//runtime: O(1)
    private double getNormalizedMutualInformation1() {
    	double MI = getMutualInformation();
    	//if MI is 0 because one partition has only 1 cluster then 
    	//the determinator is also 0 resulting in a division error.
    	//we interpret this situation as an NMI = 0
    	if(MI != 0)
    		return MI / Math.sqrt(entropies[0] * entropies[1]);
    	else
    		return 0;
    }
    
    //calculates the Fred & Jain Normalized Mutual Information for two partitions of the same size
  	//according to Wagner & Wagner (Comparing Clusterings; 2007)
  	//runtime: O(1)
    private double getNormalizedMutualInformation2() {
    	double MI = getMutualInformation();
    	//if MI is 0 because one partition has only 1 cluster then 
    	//the determinator is also 0 resulting in a division error.
    	//we interpret this situation as an NMI = 0
    	if(MI != 0)
    		return 2 * MI / (entropies[0] + entropies[1]);
    	else
    		return 0;
    }
    
    //calculates the Variation of Information for two partitions of the same size
  	//according to Gosh & Acharya (Cluster ensembles; 2011)
  	//runtime: O(1)
    private double getVariationOfInformation() {
    	return entropies[0] + entropies[1] - 2 * getMutualInformation();
    }
    
    //calculates the entropy of the partition 0 or 1 according to:
    //Gosh & Acharya (Cluster ensembles; 2011)
    //runtime: O(C)
    private double getEntropy(int partitionIndex) {
    	double entropy = 0;
    	for(int i = 0; i < clusterSizes[partitionIndex].length;i++)
    		entropy += getPLogP(clusterSizes[partitionIndex][i]);
    	return -1 * entropy;	
    }
    
    //calculates the joint entropy of the partitions 0 and 1 according to:
    //Gosh & Acharya (Cluster ensembles; 2011)
    //runtime: O(C^2)
    private double getJointEntropy() {
    	double jointEntropy = 0;
    	for(int i = 0; i < ISM.length;i++)
    		for(int j = 0; j < ISM[0].length;j++)
    			jointEntropy += getPLogP(ISM[i][j]);
    	return -1 * jointEntropy;
    }
    
    //runtime: O(1)
    private void updateEntropy(int partitionIndex, int j, int k) {
    	entropies[partitionIndex] = entropies[partitionIndex] 
		    						+ getPLogP(clusterSizes[partitionIndex][j])
		    						- getPLogP(clusterSizes[partitionIndex][j] - 1)
		    						+ getPLogP(clusterSizes[partitionIndex][k-1])
		    						- getPLogP(clusterSizes[partitionIndex][k-1] + 1);
    }
    
    //runtime: O(1)
    private void updateJointEntropy(int i, int j, int k) {
    	jointEntropy = jointEntropy 
						+ getPLogP(ISM[i][j])
						- getPLogP(ISM[i][j] - 1)
						+ getPLogP(ISM[i][k-1])
						- getPLogP(ISM[i][k-1] + 1);
    }
    
    //runtime: O(1)
    private double getPLogP(long count) {
    	return count != 0 ? (count / (double) partitions[0].length) * Math.log(count / (double) partitions[0].length) : 0;
    }
    
    //runtime: O(C^2)
    private void calcEntropies() {
    	double[] entropies = {getEntropy(0),getEntropy(1)};
    	this.entropies = entropies;
    	jointEntropy = getJointEntropy();
    }
    
	//calculates a mismatch matrix according to
	//Hubert & Arabie (Comparing Partitions; 1985)
	//runtime: O(N + C^2)
	private void calcMismatchMatrix() {		
		MM = new long[4];		
		for(int i=0;i< ISM.length;i++) {
			for(int j=0;j< ISM[0].length;j++) {
				long square = ISM[i][j] * ISM[i][j];
				MM[0] += ISM[i][j] * (ISM[i][j] - 1);
				MM[1] += square;
				MM[2] -= square;
				MM[3] -= square;
			}
			long square2 = clusterSizes[0][i] * clusterSizes[0][i];
			MM[1] -= square2;
			MM[2] += square2;
		}
		for(int j=0;j < ISM[0].length;j++) {
			long square3 = clusterSizes[1][j] * clusterSizes[1][j];
			MM[1] -= square3;
			MM[3] += square3; 
		}
		
		MM[0] /= 2; // same | same
		MM[1] = (MM[1] + partitions[0].length * partitions[0].length) / 2; // different | different
		MM[2] /= 2; // same | different
		MM[3] /= 2; // different | same
	}
	
	//calculates the size of the intersection of every cluster pair
	//runtime: O(N)
    private void calcIntersectionSizeMatrix() {
    	numOfClusters = new int[2];
    	int[][] clusterIndices = new int[2][partitions[0].length];
    	
    	//count number of clusters (and sizes) in both partitions and assign new cluster numbers starting from 1
    	for(int i=0;i<2;i++)
			for(int x=0;x<partitions[0].length;x++) {
				if (clusterIndices[i][partitions[i][x]] == 0)
					clusterIndices[i][partitions[i][x]] = ++numOfClusters[i];
				partitions[i][x] = clusterIndices[i][partitions[i][x]];
			}
		
		//set the maximum expected number of clusters. Keep this low to avoid heap size errors
		maxExpectedNumOfClusters = Math.min(Math.max(numOfClusters[0], numOfClusters[1] + 1) * 2, partitions[1].length + 1);
    	
    	clusterSizes = new int[2][maxExpectedNumOfClusters];
		
		//count the number of shared cluster members
		ISM = new long[numOfClusters[0]][maxExpectedNumOfClusters];
		for(int i = 0;i < partitions[0].length;i++) {
			ISM[partitions[0][i] - 1][partitions[1][i] - 1]++;
			clusterSizes[0][partitions[0][i] - 1]++;
			clusterSizes[1][partitions[1][i] - 1]++;
		}
    }
}
