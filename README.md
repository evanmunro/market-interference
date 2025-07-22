## Overview

This repository includes code and instructions for replicating the figures and tables
in the paper "Treatment Effects in Market Equilibrium", authored by Evan Munro, Stefan
Wager and Kuang Xu.  One main file runs all of the code to generate the data for the 3 figures and 3 tables in the paper. The replicator should expect the code to run for less than one hour.

The code requires a current [Julia](https://docs.julialang.org/en/v1/) and R installation.

## Data Availability and Provenance Statements

### Statement about Rights

- [x] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

### Summary of Availability

- [x] All data **are** publicly available.

### Details on the Data

The only data source used is the replication package for the Filmer et al. (2023). The replication package is available on the [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SGJDLC) and is in the public domain.

Filmer, Deon; Friedman, Jed; Kandpal, Eeshani; Onishi, Junko, 2021, "Replication data for: Cash Transfers, Food Prices, and Nutrition Impacts on Ineligible Children", https://doi.org/10.7910/DVN/SGJDLC, Harvard Dataverse, V1

The following files are used from this dataset:
   - balance_data_final.dta
   - moduleA-all_1.dta
   - moduleD.dta

Since we do not provide a copy of this data ourselves, download the replication package from the [Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SGJDLC)
and save these three files in the folder `data/philippines`.

## Computational requirements

### Software Requirements

The code was last run with the following versions. 

- Julia 1.10.4 
  - `Revise` 3.6.2
  - `StatFiles` 0.8.0
  - `StatsPlots` 0.15.7
  - `Optim` 1.9.4
  - `Distributions` 0.25.112
  - `Statistics` 1.10.0
  - `GLM` 1.9.0
  - `DataFrames` 1.3.6
  - `Plots` 1.40.8
  - `LaTexTabulars` 1.0.0
  - `LaTeXStrings` 1.4.0
  - `FixedEffectsModels` 1.11.0
  - `RCall` 0.14.6
  - `Roots` 2.2.5 
- R 4.4.2
  - `grf` 2.3.2

Note that as Julia releases new versions of the language and packages, and in some cases abandons old ones, it may not always be straightforward to install exactly these versions of the packages. In that case, the replicator should install the most recent version, and expect minor differences in the output. This can be accomplished by Pkg.add("PkgName"), rather than using Pkg.add(PackageSpec(name="PkgName", version = "PkgVersion")), as currently in the optional file `install_dependencies.jl`, which is described in more detail below. 

### Controlled Randomness
- Random seed for Julia is set before the generation of each table or figure in the file `replication.jl`
- Random seed for the calling of R is set on line 68 in the file `hte_analysis.jl`

Due to quirks with the `grf` package in R across operating systems, see details [here](https://grf-labs.github.io/grf/REFERENCE.html#forests-predict-different-values-depending-on-the-platform-even-though-the-seed-is-the-same), a replicator that does not use the same operating system and chip architecture should expect minor differences in Figure 2, and Tables 1-2, despite the fact that seeds are set. Specifically, Figure 2 may have an optimal treatment rule that is slightly steeper or shallower than the figure in the publication, and the distributions of points may be slightly more or less clustered. The numbers in Table 1 and Table 2 may differ by a handful of digits after the decimal point. 

### Memory, Runtime, Storage Requirements

The code needs less than one hour to run on a standard 2024 desktop or laptop machine. The storage space needed is 25 MB - 250 MB.

The code was last run on an M3 Pro with MacOS version 15.3.1 with 890GB of free space. 

### License

The code is licensed under CC0 1.0 Universal license. 

## Instructions to Replicators

#### 1. Download Data
Download the replication package of Filmer et. al (2023), available [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SGJDLC) and place the files `moduleA-all_1.dta`, `moduleD.dta`, and `balance_data_final.dta` in the `data/philippines` folder.

#### 2. Reproduce All Figures and Tables

Next, set your working directory to the `code/` folder.

If you do not already have the required Julia packages installed, first run 
``` 
julia install_dependencies.jl 
``` 

`grf` can be installed by running `install.packages("grf")` in any R console window. 

Once you have the required packages installed, run 
```
julia replicate.jl
```

The tables and figures are saved in the `exhibits/` folder.

## List of tables and figures
The script `replication.jl` produces all tables and figures in the paper using routines called from
`php_household_model.jl`.

- Figure 1 is saved as `exhibits/fig1a.pdf` and `exhibits/fig1b.pdf`
- Figure 2 is saved as `exhibits/fig2.pdf`
- Table 1 is saved as `exhibits/table1.tex`
- Table 2 is saved as `exhibits/table2.tex`
- Table 3 in the Appendix is saved as `exhibits/table3.tex`

## References

Filmer, Deon; Friedman, Jed; Kandpal, Eeshani; Onishi, Junko, 2021, "Replication data for: Cash Transfers, Food Prices, and Nutrition Impacts on Ineligible Children", https://doi.org/10.7910/DVN/SGJDLC, Harvard Dataverse, V1

Deon Filmer, Jed Friedman, Eeshani Kandpal, Junko Onishi; Cash Transfers, Food Prices, and Nutrition Impacts on Ineligible Children. The Review of Economics and Statistics 2023; 105 (2): 327–343. doi: https://doi.org/10.1162/rest_a_01061
