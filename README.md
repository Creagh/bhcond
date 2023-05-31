# Bayesian Hierarchical Clustering of Network Data (bhcond)

## Summary

`bhcond` is a [Blang](https://www.stat.ubc.ca/~bouchard/blang/) library for performing hierarchical clustering on network data. Current features focus primarily on Bayesian inference of a dendrogram from network data. This library is intended for network data that can be represented as a simple graph (i.e., an undirected, unweighted graph, with no self-loops). The dendrogram is treated as a latent parameter, representing a hierarchical clustering of the network nodes/vertices.

Four distinct model variants are available. This includes three fully Bayesian models from the class characterized by the Hierarchical Block Sampling construction [1]. This also includes the Hierarchical Random Graph model by Clauset et al. [2].  

## Installation

For installation of Blang, please refer to the [project webpage](https://www.stat.ubc.ca/~bouchard/blang/index.html).

Installing the `bhcond` library:

1. Clone this repository
2. From the command line, change to the directory with `cd bhcond`
3. Build the library using `.gradlew/build`

### Required Software

Currently, the recommended version of **Java** is **1.8.0_192**.

For producing post-processing plots, **R** version **4.0.0** or newer is recommended.

## Example

As a working example, the following code can be run from the command line. This code will fit the *Geometric-Binomial* model variant to Zachary's Karate network data [3].

```
java -cp build/install/bhcond/lib/\* bhcond.BHCOND 
--model.graph.file data/karate/karate.csv
--model.variant GEOBINOM
--engine PT 
--engine.nChains 3 
--engine.nScans 1000 
--engine.initialization COPIES 
--engine.usePriorSamples true 
--postProcessor DefaultPostProcessor 
--experimentConfigs.recordGitInfo false 
--treatNaNAsNegativeInfinity true 
--checkIsDAG false 
```

The first two Blang arguments, preceded by the command `--model`, specify the source of the data and the model variant.

The command `--engine PT` selects the non-reversible Parallel Tempering algorithm for approximate inference of the dendrogram. This is the recommended inference engine for `bhcond`.

The commands `--engine.nChains` and `--engine.nScans` control the number of annealed Markov chains and number of posterior samples to draw, respectively. Provided sufficient compute resources, the corresponding values should be increased to produce more accurate posterior approximations. Please refer to Chapter 6 and Table 6.2 from [1] for guidelines.

The following commands can also be included control the amount of thinning, and the amount of local exploration between swap attempts during the Parallel Tempering algorithm.

```
--engine.thinning 10 
--engine.nPassesPerScan 9 
```

Further details on the code arguments can be found in Section 5.1.2 of [1].

### Specifying Priors

The Blang file [`BHCOND.bl`](src/main/java/bhcond/BHCOND.bl), and specifically the `laws` block define the priors and likelihood for the model.

The code in lines 25 - 30 specify the priors and their parameters. They can be modified as desired.

https://github.com/Creagh/bhcond/blob/0285da99ee186bfa8d9e3b62b601b50ee3a4b4e9/src/main/java/bhcond/BHCOND.bl#L25-L30

However, the uniform prior over dendrograms (line 23) and likelihood (line 33) should **not** be modified.

## References

[1] Briercliffe, C. Bayesian Models for Hierarchical Clustering of Network Data. University of British Columbia. (Doctoral Dissertation), 2023.
 
[2] A. Clauset, C. Moore, and M. E. Newman. [Hierarchical structure and the prediction of missing links in networks](https://arxiv.org/pdf/0811.0484.pdf). Nature, 453(7191):98, 2008.

[3] W. W. Zachary. An information flow model for conflict and fission in small groups. Journal of anthropological research, 33(4):452–473, 1977.
