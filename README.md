## Replication Code for Treatment Effects in Market Equilibrium

This repository includes code and instructions for replicating the simulations in
[this paper](https://arxiv.org/pdf/2109.11647.pdf).

The code requires a current Julia installation and the following packages.

```
Revise, StatsPlots, Optim, LinearAlgebra, Random, Distributions, Suppressor, Statistics
```

Reproducing the simulations should take less than an hour.

#### 1. Reproducing Figure 1 and Table 1

Run the following script:
```
julia multiple_good_figure.jl
```

The components making up the table are printed. The figure is saved as `simulation.pdf`

#### 2. Reproducing Table 2

Run the following script to print the coverage results.
```
julia coverage_simulation.jl
```
This script also plots some figures comparing the estimated standard deviation to the
actual standard deviation over the monte carlo simulation runs.
