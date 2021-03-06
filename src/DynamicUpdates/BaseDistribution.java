package projects.DynamicUpdates;

/**
 * Abstract class representing the basic interface and bookkeeping mechanisms for managing
 * a base distribution of a nonparametric Bayesian model.
 * 
 * @author Thomas J. Glassen
 */

import java.util.Arrays;

public abstract class BaseDistribution {	
	
	private int[] partition;
	private int[] clusterSizes;
	private int numOfClusters;
	private int firstEmptyClusterIndex;
	
	protected void initPartition(int numOfObservations){
		partition = new int[numOfObservations];
		clusterSizes = new int[numOfObservations+1];
		Arrays.fill(partition, -1);
		numOfClusters = 0;
		firstEmptyClusterIndex = 0;
	}
	
	public void setClusterIndexOfObservation(int clusterIndex, int observationIndex) {
		if(partition[observationIndex] != -1) {
			if(clusterSizes[partition[observationIndex]]-- == 1) {
				numOfClusters--;
				if(partition[observationIndex] < firstEmptyClusterIndex)
					firstEmptyClusterIndex = partition[observationIndex];
			}
		}
		if(clusterIndex != -1) {
			if(clusterSizes[clusterIndex]++ == 0){
				numOfClusters++;
				if(clusterIndex == firstEmptyClusterIndex)
					saveFirstEmptyClusterIndex(clusterIndex + 1);
			}
		}
		partition[observationIndex] = clusterIndex;
	}
	
	private int[] getNSpecificIndices(int[] arr, int N, int conditionValue, boolean parity) {
		int[] indices = new int[N];
		int counter = 0;
		for(int i=0;i<arr.length;i++) {
			if(parity && arr[i] == conditionValue)
					indices[counter++] = i;
			else if(!parity && arr[i] != conditionValue)
					indices[counter++] = i;
			
			if(counter == N)
				break;
		}
		return indices;
	}
	
	private void saveFirstEmptyClusterIndex(int startFrom) {
		int i = startFrom;
		while(clusterSizes[i] > 0)
			i++;
		firstEmptyClusterIndex = i;
	}
	
	public int getFirstEmptyClusterIndex(){
		return firstEmptyClusterIndex;
	}
	
	protected int[] getAllObservationIndicesOfCluster(int clusterIndex) {
		return getNSpecificIndices(partition, clusterSizes[clusterIndex],clusterIndex,true);
	}
	
	public int[] getAllClusterIndices(){
		return getNSpecificIndices(clusterSizes, numOfClusters,0,false);	
	}
	
	public int getClusterIndexOfObservation(int observationIndex) {
		return partition[observationIndex];
	}
	
	public int[] getPartition(){
		return partition;
	}
	
	public int getClusterSize(int clusterIndex) {
		return clusterSizes[clusterIndex];
	}
	
	public int getNumOfObservations(){
		return partition.length;
	}

	public int getNumOfClusters(){
		return numOfClusters;
	}
		
	public void removeObservationFromCluster(int observationIndex){
		setNewClusterForObservation(-1,observationIndex);
	}
	
	public void addObservationToCluster(int clusterIndex, int observationIndex){		
		setNewClusterForObservation(clusterIndex,observationIndex);
	}
	
	public void createClusterForObservation(int observationIndex){
		setNewClusterForObservation(getFirstEmptyClusterIndex(),observationIndex);
	}
		
	public double getEmptyClusterLikelihood(int observationIndex){
		return getClusterLikelihoods(observationIndex, new int[] {getFirstEmptyClusterIndex()})[0];
	}

	public double getEmptyClusterLogLikelihood(int observationIndex){
		return getClusterLogLikelihoods(observationIndex, new int[] {getFirstEmptyClusterIndex()})[0];
	}
	
	protected abstract void setNewClusterForObservation(int newClusterIndex, int observationIndex);
	public abstract double[] getClusterLikelihoods(int observationIndex, int[] clusterIndices);
	public abstract double[] getClusterLogLikelihoods(int observationIndex, int[] clusterIndices);
}
