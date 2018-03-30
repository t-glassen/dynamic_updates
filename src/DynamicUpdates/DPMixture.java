package projects.DynamicUpdates;

/**
 * A Dirichlet Process mixture model described by a 
 * base distribution H and a concentration parameter alpha.
 * 
 * @author Thomas J. Glassen
 */

import java.util.SplittableRandom;

public class DPMixture extends MixtureModel{
	
	private double alpha, logAlpha;
	private SplittableRandom r;
	
	private boolean useLogLikelihood;
	
	public DPMixture(BaseDistribution H, double alpha){
		
		assert alpha > 0;
		
		this.H = H;
		this.alpha = alpha;
		logAlpha = Math.log(alpha);
		
		long seed = System.currentTimeMillis();
		System.out.println("[DP] seed is: " + seed);		
		r = new SplittableRandom(seed);
		
		//toggle here if log or normal likelihood should be used
		useLogLikelihood = true;
	}

	protected void iterate(){ 
		for(int obs=0;obs<H.getNumOfObservations();obs++){
			//System.out.println("current obs in iterate: " + obs);
			
			// First, remove the observation from its cluster ...
			H.removeObservationFromCluster(obs);
			int[] clusterIndices = H.getAllClusterIndices();
			
			// if there was no cluster at all, create a new one with the observation and continue with the next
			if (clusterIndices.length == 0){
				H.createClusterForObservation(obs);
				continue;
			}
			int k = 0;
			if(useLogLikelihood) {
				double[] logProbabilities = new double[H.getNumOfClusters() + 1];
				double[] clusterLogLikelihoods = H.getClusterLogLikelihoods(obs, clusterIndices);
				
				//...then get all cluster probabilities ...
				for(int j=0;j<clusterIndices.length;j++)
					logProbabilities[j] = clusterLogLikelihoods[j] + Math.log(H.getClusterSize(clusterIndices[j]));
					
				logProbabilities[clusterIndices.length] = H.getEmptyClusterLogLikelihood(obs) + logAlpha;			
				
				//... choose a new cluster ...
				k = LogExp.randCategory(logProbabilities, clusterIndices.length + 1, r);	
			}else {
				double[] probabilities = new double[H.getNumOfClusters() + 1];
				double[] clusterLikelihoods = H.getClusterLikelihoods(obs, clusterIndices);
				
				//...then get all cluster probabilities ...
				double probabilitySum = 0;
				for(int j=0;j<clusterIndices.length;j++) {
					probabilities[j] = clusterLikelihoods[j] * H.getClusterSize(clusterIndices[j]);
					probabilitySum += probabilities[j];
				}
				probabilities[clusterIndices.length] = H.getEmptyClusterLikelihood(obs) * alpha;
				probabilitySum += probabilities[clusterIndices.length];			
				
				//... choose a new cluster ...
				double treshold = r.nextDouble() * probabilitySum;
				
				//... and get the index of it ...
				double currentProbability = probabilities[0];
				while(currentProbability < treshold)
					currentProbability += probabilities[++k];
			}
			
			//...then assign the new cluster
			if (k < clusterIndices.length)
				H.addObservationToCluster(clusterIndices[k], obs);				
			else
				H.createClusterForObservation(obs);	
		}
	}
}
