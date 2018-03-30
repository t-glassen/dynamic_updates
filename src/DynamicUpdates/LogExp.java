package projects.DynamicUpdates;

/**
 * Various tool functions for logarithms.
 * 
 * @author Thomas J. Glassen
 */

import java.util.SplittableRandom;

public class LogExp {
	/*
	 * Approximation of log-gamma of a value x by Robert Sedgewick and Kevin Wayne (2011)
	 * published at http://introcs.cs.princeton.edu/java/91float/Gamma.java.html
	 */
	public static double logGamma(double x) {
		double tmp = (x - 0.5) * log(x + 4.5) - (x + 4.5);
		double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
	                     + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
	                     +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
		return tmp + log(ser * Math.sqrt(2 * Math.PI));
	}
	
	/*
	 * Fast approximation of e^x by Ankerl (2012) and developed by 
	 * Schraudolph (1999; A Fast, Compact Approximation of the Exponential Function):
	 * https://martin.ankerl.com/2007/10/04/optimized-pow-approximation-for-java-and-c-c/
	 */
	public static double expApprox(double val) {
	    final long tmp = (long) (1512775 * val + (1072693248 - 60801));
	    return Double.longBitsToDouble(tmp << 32);
	}
	
	/*
	 * Fast approximation of the natural logarithm by Ankerl (2012) based on
	 * fast exp from Schraudolph (1999; A Fast, Compact Approximation of the Exponential Function):
	 * https://martin.ankerl.com/2007/10/04/optimized-pow-approximation-for-java-and-c-c/
	 */
	public static double logApprox(double val) {
	    final double x = (Double.doubleToLongBits(val) >> 32);
	    return (x - 1072632447) / 1512775;
	}
	
	public static double exp(double val) {
		return Math.exp(val);
	}
	
	public static double log(double val) {
		return Math.log(val);
	}
	
	/*
	 * Every second cell of the array logProbabilities gives the number of repetitions of the probability
	 * in the preceding cell.
	 * We sample a category using the Gumbel-Max trick of Adams(2013) according to:
	 * https://en.wikipedia.org/wiki/Categorical_distribution#Sampling_via_the_Gumbel_distribution
	 * https://hips.seas.harvard.edu/blog/2013/04/06/the-gumbel-max-trick-for-discrete-distributions/
	 */
	public static int randCategoryFromMultipleUniforms(double[] logProbabilities, SplittableRandom r) {		
		int maxI = 0, maxC = 0, curC = 0;
		double maxGamma = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < logProbabilities.length - 1;i += 2) {
			if(logProbabilities[i + 1] == 0) //a repetition of zero marks the end of the log-probability list
				break;
			double g = -1 * log(-1 * log(r.nextDouble())); //a draw from the Gumbel distribution
			double gamma = logProbabilities[i] + log(logProbabilities[i+1]) + g;
			if(gamma > maxGamma) {
				maxGamma = gamma;
				maxI = i;
				maxC = curC;
			}
			curC += (int) logProbabilities[i+1];
		}		
		if(logProbabilities[maxI + 1] != 1) {
			int repetitions = (int) logProbabilities[maxI + 1];
			return maxC + r.nextInt(repetitions);
		}else
			return maxC;
	}
	
	/*
	 * Draws a category from a Categorical Distribution, when only log-probabilities are given
	 * using the Gumbel-Max Trick of Adams (2013):
	 * https://en.wikipedia.org/wiki/Categorical_distribution#Sampling_via_the_Gumbel_distribution
	 * https://hips.seas.harvard.edu/blog/2013/04/06/the-gumbel-max-trick-for-discrete-distributions/
	 */
	public static int randCategory(double[] logValues, int numOfValues, SplittableRandom r) {
		int c = 0;
		double maxGamma = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < numOfValues;i++) {
			double g = -1 * log(-1 * log(r.nextDouble())); //a draw from the Gumbel distribution
			double gamma = logValues[i] + g;
			if(gamma > maxGamma) {
				maxGamma = gamma;
				c = i;
			}			
		}		
		return c;
	}
	
	/*
	 * Applies the log-sum-exp trick to avoid underflows if we need 
	 * log(x_1 + x_2 + ... + x_n) while having only log(x_1) to log(x_n):
	 * https://en.wikipedia.org/wiki/LogSumExp
	 */
	public static double applyLogSumExpTrickAgainstUnderflow(double[] logValues){
		double minLogValue = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < logValues.length; i++)
			if(logValues[i] < minLogValue)
				minLogValue = logValues[i];
		
		double sumExp = 0;
		for(int i = 0; i < logValues.length; i++)
			sumExp += exp(logValues[i] - minLogValue);
		
		return log(sumExp) + minLogValue;
	}
	
	/*
	 * Calculates log(a + b), if only log(a) and log(b) are available and avoids underflows based on:
	 * https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation.2Fsubtraction
	 */
	public static double log1p(double logA, double logB) {
		if(logA > logB)
			return logB + log(1 + exp(logA - logB));
		else
			return logA + log(1 + exp(logB - logA));		
	}
}
