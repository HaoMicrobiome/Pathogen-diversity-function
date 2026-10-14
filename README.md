# Pathogen-diversity-function relationships

R scripts for three analyses of microbial diversity-function relationships along plant-pathogenic fungal abundance gradients.

## Analysis scripts

| Script | Purpose |
| --- | --- |
| `sliding_window_analysis.R` | Fits linear mixed-effects models within moving windows along the pathogen-abundance gradient to estimate changes in the microbial multidiversity-function relationship. Models include random intercepts for sampling site and tree species. |
| `mediation.R` | Estimates indirect and direct associations between pathogen abundance and diversity effects through root-effect strength, using the `mediation` package and 5,000 bootstrap replicates. |
| `lavaan_SEM_Fig3.R` | Fits structural equation models separately for low- and high-pathogen groups using `lavaan`, and exports standardized path coefficients and model-fit statistics. |

Root-effect strength describes microbial community dissimilarity between paired rhizosphere and bulk soils, averaged across bacterial and fungal communities. Diversity effects are standardized slopes relating microbial multidiversity to nutrient-acquisition functions.

## How to run

Study data are not included. Prepare the input tables described in [data/README.md](data/README.md) and place them in `data/`. Each analysis uses a separate prepared input; the sliding-window script does not generate the mediation or SEM inputs.

Install the required R packages once:

```r
install.packages(c("lme4", "lmerTest", "mediation", "lavaan", "openxlsx"))
```

Open the desired script in RStudio and click **Source**. Alternatively, with the repository as your working directory, run one analysis, for example:

```r
source("mediation.R")
```

To run all three analyses after supplying all required inputs:

```r
source("00_run_all.R")
```

Keep `analysis_helpers.R` alongside the scripts. Results are saved in timestamped folders under `outputs/` as CSV tables, with diagnostic records and R session information.
