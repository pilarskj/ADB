# ADBP

The ADBP project provides a [BEAST 2](http://www.beast2.org) package for performing phylodynamic inference under the Age-Dependent Branching Process.

## Model

Consider a binary branching process initialised with a single particle (i.e. cell) at $`t = t_{origin}`$. Each cell lives for time $`t_l`$, where the lifetimes follow a distribution with probability density function $`f`$ and cumulative density function $`F`$. At the end of its lifetime, a cell dies with probability $`d`$, or divides and generates two daughter cells. At present ($`t=0`$) , each cell is sampled with probability $`\rho`$.

Here, we consider Erlang distributed (Gamma distributed with scale $`\alpha \in \mathbb{R}_{>0}`$ and shape $`\beta \in \mathbb{N}_{>0}`$) lifetimes. Biologically, this could represent cells having $`\beta`$ transitory life stages, each with i.i.d. exponential holding times. If $`\beta = 1`$, the cell lifetimes are exponentially distributed with rate $`\alpha`$ and the process becomes a constant rate birth-death process with birth rate $`\alpha(1-d)`$ and death rate $`\alpha d`$. As $`\beta \to \infty`$, the cell divisions and deaths become synchronous and occur at regular time intervals.

Jones[^readme-1] derived the likelihood of a reconstructed phylogenetic tree $`\mathcal{T}`$ from a sample of cells under ADBP given a set of parameters $`\Theta=\{\alpha, \beta, d, \rho\}`$:

[^readme-1]: Jones G. Calculations for multi-type age-dependent binary branching processes. Journal of Mathematical Biology. 2011 Jul;63(1):33-56.

```math
\mathcal{L} (\mathcal{T} | \Theta) = \prod_{e \in \mathcal{E}} P_1(t_e) \prod_{e \in \mathcal{I}} b(s_e, t_e)
```

where $`\mathcal{E}`$ denote the external and $`\mathcal{I}`$ the internal edges (or branches) of the tree, and $`(s_e, t_e)`$ denote the start and end time of an edge $`e`$ (backwards in time). The calculation is based on integral equations $`P_0(t)`$ (probability of no sampled descendants), $`P_1(t)`$ (probability of exactly one sampled descendant) and $`b(s,t)`$ (probability of exactly one descendant at time $`s`$) for a cell originating at time $`t`$. The three equations can be solved numerically by an iterative procedure. Integrals in the equations are convolutions that can be evaluated using FFT algorithm (this involves discretizing the time intervals into equal steps).

An alternative to solving the integral equations $`b(s_e,t_e)`$ for each internal edge $`e`$ is approximating the quantity by an explicit formula (details in attached documents). This makes the likelihood calculation scalable to large trees and applicable in a Bayesian inference framework.

## Repository

The class *adbp.GammaLogLikelihood* contains all functions for solving and approximating the integral equations $`P_0(\cdot)`$, $`P_1(\cdot)`$ and $`b(\cdot,\cdot)`$, and computing the log likelihood of a phylogenetic tree based on a set of model parameters and branching times (**to review**).

The class *adbp.GammaBranchingModel* extends the BEAST 2 class *beast.base.evolution.speciation.SpeciesTreeDistribution*. Thus, ADBP can serve as a phylodynamic model (tree prior) in a MCMC chain for Bayesian phylogenetic inference (**to review**).

All the remaining classes (*adbp.Branch*, *adbp.BranchList*, *adbp.MTLogLikelihood*) attempt to extend the model to the multi-type case (**can be ignored for now**).

Also, this repository contains functions for simulating phylogenies under ADBP in R (for a fixed origin, or for a predefined number of sampled cells). An attempt to simulate such phylogenies in BEAST 2 is provided in the class *adbp.Simulator*.
