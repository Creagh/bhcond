package bhcond;

public enum ModelVariant {
	BINOM, // Beta prior on internal node probabilities and Binomial likelihood
	ZEROINFLATED, // mixture of point mass at 0 and BINOM model 
	GEOBINOM, // mixture of Truncated Geometric and BINOM model
	HRG // Clauset et al. Hierarchical Random Graph model, using maximum-likelihood estimators for internal node probs.
}
