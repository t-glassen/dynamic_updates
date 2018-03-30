package projects.DynamicUpdates;

/**
 * Provides various methods for partition distributions
 * and consensus clustering.
 * 
 * @author Thomas J. Glassen
 */

import java.util.SplittableRandom;

public class PartitionDistribution {
	
	private int[][] partitions;
	private int[] initialMP;
	private SplittableRandom r;
    
	// constructors
	/////////////////////////////////////////////////////////////////////
	
	//create an empty instance
	public PartitionDistribution() {
		partitions = null;
		initialMP = null;
		initRNG();
	}

	public PartitionDistribution(int[][] partitions){
		setNewPartitions(partitions);
		initRNG();
	}
	
	//generates K random partitions of size N with C clusters
	public PartitionDistribution(int K, int N, int C) {
		initRNG();
		setNewPartitions(K,N,C);
	}
	
	// initialization routines
	/////////////////////////////////////////////////////////////////////
	
	private void initRNG() {
		long seed = System.currentTimeMillis();
		System.out.println("[PartitionDistribution] seed is: " + seed);		
		r = new SplittableRandom(seed);
	}
	
	public void setNewPartitions(int[][] partitions) {
		this.partitions = new int[partitions.length][partitions[0].length];
		for(int i = 0; i < partitions.length;i++)
			this.partitions[i] = partitions[i].clone();
		
		initialMP = partitions[0].clone();
	}
	
	//generates K random partitions of size N with C clusters
	public void setNewPartitions(int K, int N, int C) {
		partitions = new int[K][N];		
		for(int i = 0; i < K; i++)
			partitions[i] = generateRandomPartition(N,C);
		
		initialMP = partitions[0].clone();
	}

	//sets the starting point for the local search procedures
	//(the initial mean partition)
	public void setInitialMP(int[] partition) {
		initialMP = partition.clone();
	}
	
	// random partition generation
	/////////////////////////////////////////////////////////////////////
	
	//generates a random partition of length N with C clusters
	public int[] generateRandomPartition(int N, int C) {
		int[] partition = new int[N];
		int numOfClusters = 0;
		while(numOfClusters != C) {
			numOfClusters = 0;
			int[] clusterSizes = new int[C];
			for(int i = 0; i < N; i++) {
				partition[i] = r.nextInt(C);
				if(++clusterSizes[partition[i]] == 1)
					numOfClusters++;
			}
		}
		return partition;
	}	
	
	// routines for gathering local search procedure runtime information
	/////////////////////////////////////////////////////////////////////
	
	private int[] numOfIterations;
	private int currentRun;
	
	public void prepareArrayToKeepIterations(int numOfRuns) {
		numOfIterations = new int[numOfRuns];
		currentRun = 0;
	}
	
	public int[] getNumOfIterations() {
		return numOfIterations;
	}
		
	// routines for calculating sums or averages of similarity to a given partition
	/////////////////////////////////////////////////////////////////////
	
	/*
	 * returns the sum of the specified similarity of every available partition to p
	 * runtime (transfer distance): O(K * (N + C^3))
	 * runtime (counting pairs measures): O(K * (N + C^2))
	 */
	public double getSumOfSimilarities(int[] p, int measure){
		double sumOfSimilarities = 0;
		for(int i=0;i<partitions.length;i++)
			sumOfSimilarities += getSimilarity(p,partitions[i], measure);
		
		return sumOfSimilarities;
	}
	
	/*
	 * returns the average similarity to p
	 * runtime (transfer distance): O(K * (N + C^3))
	 * runtime (counting pairs measures): O(K * (N + C^2)) 
	 */
	public double getAverageSimilarity(int[] p, int measure) {
		return getSumOfSimilarities(p, measure) / (double) partitions.length;
	}
	
	/*
	 * returns a squared matrix with the partition distances between every pair 
 	 * of the available partitions
     * runtime: O(K^2 * (N + C^3))
     */
	private int[][] getPartitionDistanceMatrix(){	
		int[][] partitionDistances = new int[partitions.length][partitions.length];
		for(int i=0;i < partitions.length;i++)			
			for(int j=i+1;j < partitions.length;j++){
				int tmpPartitionDistance = getPartitionDistance(partitions[i],partitions[j]);
				partitionDistances[i][j] += tmpPartitionDistance;
				partitionDistances[j][i] += tmpPartitionDistance;
			}		
		return partitionDistances;
	}
	
	/*
	 * returns a squared matrix with the partition similarities between every pair 
 	 * of the available partitions
 	 * runtime (transfer distance): O(K^2 * (N + C^3))
	 * runtime (counting pairs measures): O(K^2 * (N + C^2)) 
     */
	private double[][] getSimilarityMatrix(int measure){	
		double[][] partitionSimilarities = new double[partitions.length][partitions.length];
		for(int i=0;i < partitions.length;i++)			
			for(int j=i+1;j < partitions.length;j++){
				double tmpPartitionSimilarity = getSimilarity(partitions[i],partitions[j],measure);
				partitionSimilarities[i][j] += tmpPartitionSimilarity;
				partitionSimilarities[j][i] += tmpPartitionSimilarity;
			}		
		return partitionSimilarities;
	}
	
	/*
	 * returns the similarity between two partitions calculated via a specified measure
	 * runtime (transfer distance): O(N + C^3)
	 * runtime (counting pairs measures): O(N + C^2)
	 */
	public double getSimilarity(int[] p1, int[] p2, int measure) {
		if(measure == 0)
			return getPartitionDistance(p1,p2);
		else if(measure == 11)
			return 1 - getPartitionDistance(p1,p2) / (double) (p1.length - 1);
		else
			return new DynamicPartitionSimilarity(p1,p2).getSimilarity(measure);
	}
	
	/*
	 * returns true if the first given similarity is better than the second,
	 * when using the specified measure.
	 * runtime: O(1)
	 */
	public boolean isBetterThan(double similarity1, double similarity2, int measure) {
		return (measure == 0 || measure == 5 || measure == 13 || measure == 15) ? similarity1 < similarity2 : similarity1 > similarity2;
	}
	
	/*
	 * returns the worst similarity on the specified measure
	 * runtime: O(1)
	 */
	public double getWorstSimilarity(int measure) {
		return (measure == 0 || measure == 5 || measure == 13 || measure == 15) ? Double.MAX_VALUE : Double.MIN_VALUE;
	}
	
	/*	Calculation of the squared partition distance via a reduction of the problem 
	 *  to a linear sum assignment problem according to Konovalov, Litow & Bajema
	 *  (Partition-distance via the assignment problem; 2005, S. 2464)
	 *  runtime: O(N + C^3)
	 */ 
	public int getPartitionDistance(int[] partition1, int[] partition2){
				
		int numOfClustersP1 = 0, numOfClustersP2 = 0;
		int[] clusterIndicesP1 = new int[partition1.length], clusterIndicesP2 = new int[partition2.length];
		
		//count number of clusters in both partitions and assign new cluster numbers starting from 1
		for(int x=0;x<partition1.length;x++){
			if (clusterIndicesP1[partition1[x]] == 0){
				numOfClustersP1++;
				clusterIndicesP1[partition1[x]] = numOfClustersP1;
			}
			if (clusterIndicesP2[partition2[x]] == 0){
				numOfClustersP2++;
				clusterIndicesP2[partition2[x]] = numOfClustersP2;
			}
		}
		
		//switch partitions if partiton 1 has more clusters 
		if(numOfClustersP2 < numOfClustersP1){
			int[] tmp = clusterIndicesP1;
			clusterIndicesP1 = clusterIndicesP2;
			clusterIndicesP2 = tmp;
			int tmp2 = numOfClustersP1;
			numOfClustersP1 = numOfClustersP2;
			numOfClustersP2 = tmp2;
		    tmp = partition1;
		    partition1 = partition2;
		    partition2 = tmp;
		}
		
		double[][] costMatrix = new double[numOfClustersP1][numOfClustersP2];
		int[] clusterSizes = new int[numOfClustersP1];
		
		//count the size of the clusters in partition 1 and decrement every cell in the cost matrix
		//where an item is shared between the corresponding clusters
		for(int x=0;x<partition1.length;x++){						
			int i = clusterIndicesP1[partition1[x]] - 1;			
			clusterSizes[i]++;
			costMatrix[i][clusterIndicesP2[partition2[x]] - 1]--;
		}
		
		//add the cluster size of cluster i to every cell in row i
		for(int i=0;i<numOfClustersP1;i++)
			for(int j=0;j<numOfClustersP2;j++)
				costMatrix[i][j] += clusterSizes[i];		
		
		int jobIndexForWorker[] = (new HungarianAlgorithm(costMatrix)).execute();
		
		int overallCosts = 0;
		for(int worker=0;worker < jobIndexForWorker.length;worker++)
			overallCosts += costMatrix[worker][jobIndexForWorker[worker]];
		
		return overallCosts;
	}
	
	// routines for calculating the lower similarity bound
	/////////////////////////////////////////////////////////////////////
	
	/*
	 * returns the minimum number of disagreements of a partition to
	 * the ensemble calculated via the lower bound given by:
	 * Bertolacci & Wirth (Are approximation algorithms for consensus clustering worthwhile?; 2007)
	 * Goder & Filkov (Consensus Clustering Algorithms: Comparison and Refinement; 2007)
	 * runtime: O(N^2 * K)
	 */
	public long getMinNumOfDisagreements() {
		long minDisagreements = 0;
		for(int i = 0;i < partitions[0].length;i++) {
			for(int j = 0;j < partitions[0].length;j++) {
				long same = 0;
				for(int k = 0;k < partitions.length;k++) {
					same += (partitions[k][i] == partitions[k][j] ? 1 : 0);
				}
				minDisagreements += (same > partitions.length - same ? partitions.length - same : same);
			}
		}
		return minDisagreements;
	}
	
	/*
	 * returns the maximum number of agreements of a partition to
	 * the ensemble calculated via the lower bound on the number of disagreements.
	 * runtime: O(N^2 * K)
	 */
	public long getMaxNumOfAgreements() {
		return partitions[0].length * (partitions[0].length - 1) / 2 - getMinNumOfDisagreements();
	}
	
	/*
	 * calculates the upper bound of the average rand index of a
	 * partition to the ensemble based on the maximum number of agreements.
	 * runtime: O(N^2 * K)
	 */
	public double getUpperBoundOfAverageRandIndex() {
		return getMaxNumOfAgreements() /(partitions.length * partitions[0].length * (partitions[0].length - 1) / 2);
	}
	
	/*
	 * calculates the lower bound of the average mirkin metric of a
	 * partition to the ensemble based on the minimum number of disagreements.
	 * runtime: O(N^2 * K)
	 */
	public double getLowerBoundOfAverageMirkinMetric() {
		return 2 * getMinNumOfDisagreements() / (double) partitions.length;
	}
	
	/*
	 * calculates the upper bound of the average bounded mirkin similarity of a
	 * a partition to the ensemble based on the minimum number of disagreements.
	 * runtime: O(N^2 * K)
	 */
	public double getUpperBoundOfAverageBoundedMirkinSimilarity() {
		return 1 - getLowerBoundOfAverageMirkinMetric() / (partitions[0].length * partitions[0].length);
	}
	
	// routines for calculating the mean, median and mode partitions
	/////////////////////////////////////////////////////////////////////
	
	//returns all median partitions
	//runtime: O(K^2 * (N + C^3))
	public int[][] getMedianPartitions(){
		int[][] partitionDistances = getPartitionDistanceMatrix();
		int[] selectedPartitionIndices = new int[partitions.length];
		
		int lowestSum = Integer.MAX_VALUE;
		int numOfMedianPartitions = 0;
		for(int i=0;i<partitionDistances.length;i++) {
			int tmpDistanceSum = 0;
			for(int j=0;j<partitionDistances[0].length;j++)
				tmpDistanceSum += partitionDistances[i][j];
			if(tmpDistanceSum < lowestSum) {
				lowestSum = tmpDistanceSum;
				numOfMedianPartitions = 1;
				selectedPartitionIndices[0] = i;
			}else if(tmpDistanceSum == lowestSum) {
				numOfMedianPartitions++;
				selectedPartitionIndices[numOfMedianPartitions - 1] = i;
			}
		}
		
		int[][] selectedPartitions = new int[numOfMedianPartitions][];
		for(int i=0;i < numOfMedianPartitions;i++)
			selectedPartitions[i] = partitions[selectedPartitionIndices[i]].clone();
		return selectedPartitions;
	}
	
	//returns first median partition
	//runtime: O(K^2 * (N + C^3))
	public int[] getMedianPartition() {
		return getMedianPartitions()[0];
	}
	
	//returns all mode partitions
	//runtime: O(K^2 * (N + C^3))
	public int[][] getModePartitions(int toleranceDistance){
		int[][] partitionDistances = getPartitionDistanceMatrix();
		int[] selectedPartitionIndices = new int[partitions.length];
		
		int highestNumOfRepetitions = 0;
		int numOfModePartitions = 0;
		for(int i=0;i<partitionDistances.length;i++) {
			int tmpRepetitions = 0;
			for(int j=0;j<partitionDistances[0].length;j++)
				if(partitionDistances[i][j] <= toleranceDistance)
					tmpRepetitions++;
			if(tmpRepetitions > highestNumOfRepetitions) {
				highestNumOfRepetitions = tmpRepetitions;
				numOfModePartitions = 1;
				selectedPartitionIndices[0] = i;
			}else if(tmpRepetitions == highestNumOfRepetitions) {
				numOfModePartitions++;
				selectedPartitionIndices[numOfModePartitions - 1] = i;
			}
		}
		
		int[][] selectedPartitions = new int[numOfModePartitions][];
		for(int i=0;i < numOfModePartitions;i++)
			selectedPartitions[i] = partitions[selectedPartitionIndices[i]].clone();
		return selectedPartitions;
	}
	
	//returns first mode partition
	//runtime: O(K^2 * (N + C^3))
	public int[] getModePartition(int toleranceDistance) {
		return getModePartitions(toleranceDistance)[0];
	}
	
	/*
	 * choose a the best of K partitions according to Goder & Filkov
	 * (Consensus Clustering Algorithms: Comparison and Refinement; 2007)
	 * runtime (transfer distance): O(K^2 * (N + C^3))
	 * runtime (counting pairs measures): O(K^2 * (N + C^2))
	 * Approximation factor: 2
	 */
	public int[] bestOfK(int measure) {		
		double[] similaritySums = new double[partitions.length];
		for(int i=0;i < partitions.length;i++)			
			for(int j=i+1;j < partitions.length;j++){
				double tmpPartitionSimilarity = getSimilarity(partitions[i],partitions[j],measure);
				similaritySums[i] += tmpPartitionSimilarity;
				similaritySums[j] += tmpPartitionSimilarity;
			}
		double bestSimilaritySum = similaritySums[0];
		int bestPartition = 0;
		for(int i=1;i < partitions.length; i++)
			if(isBetterThan(similaritySums[i], bestSimilaritySum, measure)) {
				bestSimilaritySum = similaritySums[i];
				bestPartition = i;
			}
		return partitions[bestPartition];	
	}
	
	/*
	 * combination of Pick-A-Cluster and CC-Pivot, which has 11/7 approximation
	 * according to Ailon et al. (Aggregating Inconsistent Information: Ranking
	 * and Clustering; 2005)
	 * runtime: O(K*N*C) for the whole algorithm
	 * Approximation factor: 11/7
	 */
	public int[] combinedPickAClusterCCPivot() {
		int[] solution1 = pickACluster();
		int[] solution2 = CCPivot();
		if(getSumOfSimilarities(solution1,10) > getSumOfSimilarities(solution2,10))
			return solution2;
		else
			return solution1.clone();
	}
	
	/*
	 * choose a partition at random according to Bertolacci & Wirth
	 * (Are approximation algorithms for consensus clustering worthwhile?; 2007)
	 * runtime: O(1)
	 * Approximation factor: 2
	 */
	public int[] pickACluster() {
		return partitions[r.nextInt(partitions.length)];
	}
	
	/*
	 * CC-Pivot acording to Ailon et al. (Aggregating Inconsistent Information: 
	 * Ranking and Clustering; 2005) and Bertolacci & Wirth
	 * (Are approximation algorithms for consensus clustering worthwhile?; 2007)
	 * runtime: O(K*N*C) for the whole algorithm
	 * Approximation factor: 2
	 */
	public int[] CCPivot() {
		
		int[] individualIndices = new int[initialMP.length];
		int[] clusterLabels = individualIndices.clone();
		int numOfClusters = 0;
		
		//store initial individual indices
		for(int i = 0; i < initialMP.length;i++)
			individualIndices[i] = i;
		
		int i = 0;
		while(i < initialMP.length) {
			//choose pivot
			int pivot = i + r.nextInt(initialMP.length - i);
			//switch positions
			int temp = individualIndices[i];
			individualIndices[i] = individualIndices[pivot];
			individualIndices[pivot] = temp;
			pivot = i;
			
			//set pivot into new cluster
			clusterLabels[pivot] = numOfClusters++;
			
			//search other cluster members
			int j = ++i;
			while(j < initialMP.length) {
				int wPlus = 0, wMinus = 0;
				for(int k = 0; k < partitions.length;k++) {
					if(partitions[k][individualIndices[pivot]] == partitions[k][individualIndices[j]])
						wPlus++;
					else
						wMinus++;					
				}
				//switch positions
				if(wPlus >= wMinus) {
					int temp2 = individualIndices[i];
					individualIndices[i] = individualIndices[j];
					individualIndices[j] = temp2;
					clusterLabels[i++] = clusterLabels[pivot];
				}
				j++;
			}
		}
		
		//build consensus
		int[] consensus = new int[initialMP.length];
		for(i = 0; i < initialMP.length;i++)
			consensus[individualIndices[i]] = clusterLabels[i];
		
		numOfIterations[currentRun++] = 1;
		
		return consensus;
	}
	
	/*
	 * Various similarity measures via the local optimization algorithm in a static fashion
	 * runtime: O(K*N*C*(N + C^2)) per cycle
	 */
	public int[] getMeanPartitionStatSim(int measure){		
		double bestSimilaritySum = 0;		
		
		DynamicPartition mp = new DynamicPartition(initialMP);
		bestSimilaritySum = getSumOfSimilarities(mp.getPartition(), measure);
	    	    
		boolean success = true;
	    while(success){
	        success = false;
	    	for(int obj=0;obj<initialMP.length;obj++){
	    		int[] clusterIndices = mp.getAscendingClusterNrs();	    		
	    		int oldObjClusterIndex = mp.getClusterIndexOfObj(obj), bestObjClusterIndex = oldObjClusterIndex;
	    		int sizeOfClusterOfObj = mp.getSizeOfClusterOfObj(obj);
	    		for(int clusterIndex=0;clusterIndex<clusterIndices.length;clusterIndex++){
	    			if (clusterIndices[clusterIndex] == oldObjClusterIndex){ //create new cluster
	    				if (sizeOfClusterOfObj == 1) //object is already in its own cluster
	    					continue;
	    				mp.movObjToNew(obj);
	    			}
	    			else //choose existing cluster
	    				mp.movObjTo(obj, clusterIndices[clusterIndex]);
	    			
	    			double similaritySum = getSumOfSimilarities(mp.getPartition(), measure);
	    				
	                if(isBetterThan(similaritySum,bestSimilaritySum,measure)){
	                	bestSimilaritySum = similaritySum;
	                	bestObjClusterIndex = mp.getClusterIndexOfObj(obj);
	                	success = true;
	                }
	    		}
	    		mp.movObjTo(obj, bestObjClusterIndex);
	    	}
	    }
	    return mp.getPartition();
	}
	
	/*
	 * Various similarity measures via the local optimization algorithm in a dynamic fashion
	 * runtime: O(K*N*C) per cycle
	 */
	public int[] getMeanPartitionDynSim(int measure){		
		
		double bestSimilaritySum = 0;
		
		DynamicPartitionSimilarity[] PS = new DynamicPartitionSimilarity[partitions.length];
		for(int i = 0; i < partitions.length;i++) {
			PS[i] = new DynamicPartitionSimilarity(partitions[i], initialMP);
			bestSimilaritySum += PS[i].getSimilarity(measure);
		}
	    	    
		boolean success = true;
		int numOfIter = 0;
	    while(success){
	    	numOfIter++;
	        success = false;
	    	for(int obj=0;obj<initialMP.length;obj++){
	    		int[] clusterIndices = PS[0].getAscendingClusterNrs();	    		
	    		int oldObjClusterNr = PS[0].getClusterNrOfObj(obj), bestObjClusterNr = oldObjClusterNr;
	    		int sizeOfClusterOfObj = PS[0].getSizeOfClusterOfObj(obj);
	    		for(int clusterIndex=0;clusterIndex<clusterIndices.length;clusterIndex++){
	    			double similaritySum = 0;			
	    			if (clusterIndices[clusterIndex] == oldObjClusterNr){ //create new cluster
	    				if (sizeOfClusterOfObj == 1) //object is already in its own cluster
	    					continue;
	    				for(int p = 0;p < partitions.length;p++)
	    					PS[p].movObjToNew(obj);
	    			}
	    			else { //choose existing cluster
	    				for(int p = 0;p < partitions.length;p++)
	    					PS[p].movObjTo(obj, clusterIndices[clusterIndex]);
	    			}
	    			
	    			for(int p = 0;p < partitions.length;p++)
	    				similaritySum += PS[p].getSimilarity(measure);
	    				
	                if(isBetterThan(similaritySum,bestSimilaritySum,measure)){
	                	bestSimilaritySum = similaritySum;
	                	bestObjClusterNr = PS[0].getClusterNrOfObj(obj);
	                	success = true;
	                }
	    		}
	    		for(int i = 0; i < partitions.length;i++)				
    				PS[i].movObjTo(obj, bestObjClusterNr);
	    	}
	    }
	    System.out.println("iter (" + measure + "): " + numOfIter);
	    numOfIterations[currentRun++] = numOfIter;
	    return PS[0].getPartition();
	}
	
	/* 
	 * The heuristic method of Huelsenbeck & Andolfatto 
	 * (Inference of Population Structure Under a Dirichlet Process Model; 2007, S.1792)
	 * with the faster reduction of Glassen
	 * (Psychologisch orientierte Kategorisierung in der kognitiven Robotik mit dem Hierarchischen Dirichlet Prozess; 2018)
	 * runtime: O(K*N^2*C + K*N*C*(N + C^3)) per cycle
	 */		
	public int[] getMeanPartition1(){		
		double bestDistanceSum;
		
		DynamicPartition mp = new DynamicPartition(initialMP);
		bestDistanceSum = getSumOfSimilarities(mp.getPartition(), 0);
		
		boolean success = true;
		int numOfIter = 0;
	    while(success){
			numOfIter++;
	        success = false;
	    	for(int obj=0;obj<initialMP.length;obj++){
	    		int[] clusterIndicesMP = mp.getAscendingClusterNrs();
	    		int bestObjClusterIndex = mp.getClusterIndexOfObj(obj), oldObjClusterIndex = bestObjClusterIndex;
	    		for(int clusterIndex=0;clusterIndex<clusterIndicesMP.length;clusterIndex++){
	    			if (clusterIndicesMP[clusterIndex] == oldObjClusterIndex){
						if (clusterIndicesMP.length == initialMP.length)
							continue;
						mp.movObjToNew(obj);
	    			}
	    			else
	    				mp.movObjTo(obj,clusterIndicesMP[clusterIndex]); 
 	                
	                double tmpDistanceSum = getSumOfSimilarities(mp.getPartition(), 0); //one can trade space for time via keeping the cost matrizes and updating them constantly  

	                if(tmpDistanceSum < bestDistanceSum){
	                	bestDistanceSum = tmpDistanceSum;
	                	bestObjClusterIndex = mp.getClusterIndexOfObj(obj);
	                	success = true;
	                }
	    		}
	    		mp.movObjTo(obj, bestObjClusterIndex);	    		
	    	}
	    }
		System.out.println("iter (0): " + numOfIter);
		numOfIterations[currentRun++] = numOfIter;
	    return mp.getPartition();
	}
}
