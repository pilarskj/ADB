# ADB

This repository provides a [BEAST 2](http://www.beast2.org) package for performing phylodynamic inference under the Age-Dependent Branching process.

The [preprint](https://doi.org/10.1101/2025.08.28.672870) describing the phylodynamic model and its applications is available on bioRxiv.

## The process

Consider a binary branching process initialised with a single particle (e.g. cell) at $`t = t_{origin}`$. Each cell lives for some time, where the lifetimes follow a distribution with probability density function $`f`$ and cumulative density function $`F`$. At the end of its lifetime, a cell dies with probability $`d`$, or divides and generates two daughter cells. At present ($`t=0`$), each extant cell is sampled with probability $`\rho`$.

Here, we consider Gamma- or Erlang-distributed lifetimes with shape $`k`$ and scale $`\theta`$ (mean lifetime $`\ell = k \theta`$). If $`k=1`$, the cell lifetimes are exponentially distributed with rate $`\theta`$ and the process becomes a constant rate birth-death process with birth rate $`\theta (1-d)`$ and death rate $`\theta d`$. As $`k \to \infty`$, the branching events become more regular and synchronous.

## Phylodynamic likelihood

The likelihood of a phylogenetic tree $`\mathcal{T}`$ under the age-dependent branching process has been initially derived by Jones[^readme-1].

This package provides an implementation of this likelihood and a scalable approximation which enables phylodynamic inference under the ADB process in the Bayesian inference framework.

[^readme-1]: Jones G. Calculations for multi-type age-dependent binary branching processes. Journal of Mathematical Biology. 2011 Jul;63(1):33-56.

## To Do:
- example (xml) / tutorial / installation guide
- extensions 
