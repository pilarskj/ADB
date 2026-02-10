# ADB

This repository provides a [BEAST 2](http://www.beast2.org) package for performing phylodynamic inference under the **Age-Dependent Branching** process :deciduous_tree: :chart_with_upwards_trend:.

The [preprint](https://doi.org/10.1101/2025.08.28.672870) describing the phylodynamic model and its applications is available on bioRxiv[^readme-1].

[^readme-1]: Mulberry, N., Pilarski, J., Dinger, J., & Stadler, T. (2025). Bayesian phylodynamics for developmental biology: incorporating age-dependence. *bioRxiv*, 2025-08.

## The process

Consider a binary branching process initialised with a single particle (e.g. cell) at $`t = t_{origin}`$. 
Each cell lives for some time, where the lifetimes follow a distribution with probability density function $`f`$ and cumulative density function $`F`$. 
At the end of its lifetime, a cell dies with probability $`d`$, or divides and generates two daughter cells. 
At present ($`t=0`$), each extant cell is sampled with probability $`\rho`$.

Here, we consider Gamma- or Erlang-distributed lifetimes with shape $`k`$ and scale $`\theta`$ (mean lifetime $`\ell = k \theta`$). 
If $`k=1`$, the cell lifetimes are exponentially distributed with rate $`\theta`$ and the process becomes a constant rate birth-death process with birth rate $`\theta (1-d)`$ and death rate $`\theta d`$. 
As $`k \to \infty`$, the branching events become more regular and synchronous.

## Phylodynamic likelihood

The likelihood of a phylogenetic tree $`\mathcal{T}`$ under the age-dependent branching process has been initially derived by Jones (2011)[^readme-2].

This package provides an implementation of this likelihood and a scalable approximation which enables phylodynamic inference under the ADB process in the Bayesian inference framework.

Note that the likelihood equation is solved numerically, and comes with significant computational costs over the analytical birth-death likelihood. 
We thus recommend to only use the ADB package when there is prior expectation of age-dependence :hourglass_flowing_sand:.

[^readme-2]: Jones, G. (2011). Calculations for multi-type age-dependent binary branching processes. *Journal of Mathematical Biology*, 63(1), 33-56.

## :open_book: User guide

### Software requirements

This package requires Java 17 and at least BEAST v2.7.

### Running on command line

Currently, you can run a bundled version of the package and its dependencies on the command line.

To run ADB with the example XML, use the following command:

```bash
java -jar bin/ADB.jar example/inference.xml
```
and specify the package versions using `-version_file version.xml`.

:female_detective: We recommend running the package with option `-loglevel debug`. This enables monitoring numerical errors in the phylodynamic likelihood calculation.

### Features

The `GammaBranchingModel` class defines the phylodynamic model. It takes as input a tree $`\mathcal{T}`$, parameters $`k, \ell, d, \rho`$, and $`t_{origin}`$.

In the XML, you can adjust settings for the numerical solution of the phylodynamic likelihood. The most important are:

- `approx`: if `true`, the approximation is used (in this case, $`k`$ must be an integer)
- `stepSizeP`: determines the resolution of FFTs (must be a power of $`2`$)
- `conditionOnRoot`: if `true`, the likelihood is conditioned on survival of the process since $`t_{root}`$, otherwise $`t_{origin}`$

:bulb: Using the approximation and reducing the step size massively speed up computations, at the cost of accuracy. 

In practice, we recommend to always use the approximation when running MCMC. 
We advice to use a step size of at least $`2^{10}`$ (`1024`) if $`\rho \gtrsim 50\%`$, $`2^{12}`$ (`4096`) if $`\rho \gtrsim 10\%`$, and $`2^{14}`$ (`16384`) if $`\rho \gtrsim 1\%`$.
At lower sampling proportions and when $`d`$ approaches $`0.5`$, the analyses tend to accumulate numerical errors and have convergence issues.
Additionally, if the shape parameter is expected to be low $`k \lesssim 20`$, we suggest putting informative priors on at least one of $`d, \rho`$ due to non-identifiability concerns.

:warning: Note that the model is currently validated for phylodynamic inference from *fixed* phylogenetic trees.
Operators enabling joint tree and parameter inference from sequence alignments are under development. 