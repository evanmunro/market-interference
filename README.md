## Replication Code for Treatment Effects in Market Equilibrium 

This repository includes code and instructions for replicating the figures and tables
in the paper "Treatment Effects in Market Equilibrium", authored by Evan Munro, Stefan
Wager and Kuang Xu.  

The code requires a current [Julia](https://docs.julialang.org/en/v1/) and R installation. The following Julia packages are required:

```
StatFiles, Revise, StatsPlots, Optim, LinearAlgebra, Random, Distributions, Statistics,
GLM, DataFrames, Plots, StatsPlots, LaTeXTabulars, LaTeXStrings, FixedEffectModels, RCall,
Roots
```

The following R packages are required:
```
grf
```

For the paper, figures were generated using Julia version 1.9.1 and R version 4.2.1. The replication of the tables and figures in the papers requires two steps.

#### 1. Download Data
Download the replication package of Filmer et. al (2023), available [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SGJDLC) and place the files `moduleA-all_1.dta`, `moduleD.dta`, and `balance_data_final.dta` in the `data/philippines` folder.

Reproducing the tables and figures using the below steps should take less than 1 hour on a standard computer.

#### 2. Reproduce All Figures and Tables

Set your working directory to the `code/` folder. From there, run
```
julia replicate.jl
```

The tables and figures are saved in the `exhibits/` folder.
