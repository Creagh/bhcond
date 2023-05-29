# Bayesian Hierarchical Clustering of Network Data (bhcond)

## Summary

`bhcond` is a [Blang](https://www.stat.ubc.ca/~bouchard/blang/) library for performing hierarchical clustering on network data. Current features focus primarily on Bayesian inference of a dendrogram from network data. This library is intended for network data that can be represented as a simple graph (i.e., an undirected, unweighted graph, with no self-loops).

Four distinct model variants are available. This includes three fully Bayesian models from the class characterized by the Hierarchical Block Sampling construction [1]. This also includes the Hierarchical Random Graph model by Clauset et al. [2].  


# References

[1] Briercliffe, C. Bayesian Models for Hierarchical Clustering of Network Data. University of British Columbia. (Doctoral Dissertation), 2023.
 
[2] A. Clauset, C. Moore, and M. E. Newman. [Hierarchical structure and the prediction of missing links in networks](https://arxiv.org/pdf/0811.0484.pdf). Nature, 453(7191):98, 2008.



