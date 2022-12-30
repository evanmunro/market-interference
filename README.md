## Replication Code for Treatment Effects in Market Equilibrium Paper

This repository includes code and instructions for replicating the figures and tables
in the paper "Treatment Effects in Market Equilibrium", authored by Evan Munro, Stefan
Wager and Kuang Xu.  

The code requires a current Julia and R installation. The following Julia packages are required:

```
Revise, StatsPlots, Optim, LinearAlgebra, Random, Distributions, Suppressor, Statistics,
FixedEffectModels, DataFrames, Plots
```

The following R packages are required:
```
readstata13, grf, pracma
```

Download the replication package of Gertler et. al (2012), available [here](https://www.openicpsr.org/openicpsr/project/116375/version/V1/view) and place the file `investments_data.dta` in the `data` folder.

Reproducing the tables and figures using the below steps should take less than 3 hours.

#### 1. Reproducing Figure 1 and Table 1

Run the following script:
```
julia generate_simulation_plots.jl
```

The components making up the table are printed. The figures are saved as pdfs.

#### 2. Reproducing Table 2 and Figure 2

To generate the data for Table 2 (printed) and some data for Figure 2, then run:
```
Rscript hte_analysis.R
```

Then, to generate Figure 2, run the following script:

```
julia  empirical_plots.jl
```
