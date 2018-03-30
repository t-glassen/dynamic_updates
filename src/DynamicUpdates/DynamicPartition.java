package projects.DynamicUpdates;

/**
 * This class maintains all information about a partition
 * and keeps them up to date if an individual is moved from
 * one cluster into another.
 * 
 * @author Thomas J. Glassen
 */

/*
 * This version has a worst case runtime of the unique cluster nr search in O(C), 
 * where C is the highest reached number of clusters in the whole optimization procedure  
 */

public class DynamicPartition {

	private int[] partition;
	private int numOfClusters;
	private int[] clusterSizes;
	private int firstEmptyCluster;
	private int maxExpectedNumOfClusters;
	
	public DynamicPartition(int[] partition){
		this.partition = partition.clone();
		setBookkeepingInformation();
		firstEmptyCluster = numOfClusters + 1;
	}
	
	//runtime: O(C)
	private void increaseSizeOfArrays() {		  
		//set the maximum expected number of clusters. Keep this low to avoid heap size errors
		maxExpectedNumOfClusters = Math.min((numOfClusters + 1) * 2, partition.length + 1);
		 
		int[] newClusterSizes = new int[maxExpectedNumOfClusters];
		
	  	for(int i = 0; i < clusterSizes.length;i++)
	  		newClusterSizes[i] = clusterSizes[i];
	  	clusterSizes = newClusterSizes;
	}
	
	//calculates the size of the intersection of every cluster pair
	//runtime: O(N)
    private void setBookkeepingInformation() {    	
    	numOfClusters = renumberPartition(partition);
    	
    	//set the maximum expected number of clusters. Keep this low to avoid heap size errors
    	maxExpectedNumOfClusters = Math.min((numOfClusters + 1) * 2, partition.length + 1);
    	
    	clusterSizes = new int[maxExpectedNumOfClusters];
		
		//calculate cluster sizes
		for(int i = 0;i < partition.length;i++)
			clusterSizes[partition[i]-1]++;
    }
    
    //runtime: O(N)
    private int renumberPartition(int[] partition) {
    	int numOfClusters = 0;
    	int[] clusterIndices = new int[partition.length];
    	
    	//count number of clusters and assign new cluster numbers starting from 1
    	for(int i=0;i<partition.length;i++) {
    		if (clusterIndices[partition[i]] == 0)
    			clusterIndices[partition[i]] = ++numOfClusters;
    		partition[i] = clusterIndices[partition[i]];
    	}
    	return numOfClusters;
    }
    
    //runtime: O(1)
	public void movObjTo(int index, int k){		
		//get relevant index
		int j = partition[index] - 1;
		  
		if(j == k - 1)
			return; //nothing to do		
		
    	//update cluster index
    	partition[index] = k;
    	
    	//update cluster sizes
		clusterSizes[j]--;
		clusterSizes[k-1]++;
		
		//update number of clusters
		if(clusterSizes[j] == 0) {
			numOfClusters--;
			//keep cluster nr if j + 1 is lower than firstEmptyCluster
			if(firstEmptyCluster > j + 1)
				firstEmptyCluster = j + 1;
		}
		if(clusterSizes[k-1] == 1)
			numOfClusters++;
	}
	
	//runtime: O(C)
	public void movObjToNew(int index){			  
		//get the next empty cluster. Runtime: O(C)
		int emptyClusterNr = firstEmptyCluster;
		firstEmptyCluster++;
		  
		while(firstEmptyCluster <= clusterSizes.length && clusterSizes[firstEmptyCluster - 1] > 0)
			firstEmptyCluster++;
		  
		//increase size of arrays if maxExpectedNumOfClusters is nearly reached
		if(firstEmptyCluster == maxExpectedNumOfClusters)
			increaseSizeOfArrays();
		  
		movObjTo(index, emptyClusterNr);
	}
	
	public int[] getAscendingClusterNrs(){
		int counter = numOfClusters;
		int[] clusterNrs = new int[counter];
		for(int i=0;counter > 0;i++)
			if(clusterSizes[i] > 0){
				clusterNrs[clusterNrs.length - counter] = i+1;
				counter--;
			}
		return clusterNrs;
	}
	
	//runtime: O(1)
	public int getClusterIndexOfObj(int indexOfObj){
		return partition[indexOfObj];
	}
	
	//runtime: O(1)
	public int getSizeOfClusterOfObj(int indexOfObj){
		return clusterSizes[partition[indexOfObj] - 1];
	}
	
	//runtime: O(1)
	public int[] getPartition(){
		return partition;
	}
	
	//runtime: O(N)
	public int[] getRenumberedPartition() {
		int[] newPartition = partition.clone();
		renumberPartition(newPartition); 
		return newPartition;
	}
}
