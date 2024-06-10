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

The citation for the data is as follows:
Filmer, Deon; Friedman, Jed; Kandpal, Eeshani; Onishi, Junko, 2021, "Replication data for: Cash Transfers, Food Prices, and Nutrition Impacts on Ineligible Children", https://doi.org/10.7910/DVN/SGJDLC, Harvard Dataverse, V1

The following files are used from this dataset:
   - balance_data_final.dta
   - moduleA-all_1.dta
   - moduleD.dta

Since we do not provide a copy of this data ourselves, download the replication package from the [Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SGJDLC)
and save these three files in the folder `data/philippines`.

## Computational requirements

### Software Requirements

- Julia 1.9.1
  - `Revise` 3.5.3
  - `StatFiles` 0.8.0
  - `StatsPlots` 0.15.5
  - `Optim` 1.7.6
  - `Distributions` 0.25.95
  - `Statistics` 1.9.0
  - `GLM` 1.8.3
  - `DataFrames` 1.5.0
  - `Plots` 1.38.15
  - `LaTexTabulars` 0.1.3
  - `LaTeXStrings` 1.3.1
  - `FixedEffectsModels` 1.11.0
  - `RCall` 0.13.15
- R 4.2.1
  - `grf` 2.2.0

### Controlled Randomness
- Random seed for Julia is set before the generation of each table or figure in the file `replication.jl`
- Random seed for the calling of R is set on line 68 in the file `hte_analysis.jl`


### Memory, Runtime, Storage Requirements

The code needs less than one hour to run on a standard 2024 desktop or laptop machine. The storage space needed is 25 MB - 250 MB.

The code was last run on a 4-core Intel-based laptop with MacOS version 13.6.7 with 47GB of free space.

### License

The code is licensed under a TBD license.

## Instructions to Replicators

#### 1. Download Data
Download the replication package of Filmer et. al (2023), available [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SGJDLC) and place the files `moduleA-all_1.dta`, `moduleD.dta`, and `balance_data_final.dta` in the `data/philippines` folder.

#### 2. Reproduce All Figures and Tables

Set your working directory to the `code/` folder. From there, run
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
