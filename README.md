# Selected statistical analyses

Selected R analyses supporting the study of microbial diversity-function relationships along putative plant-pathogenic fungal abundance gradients.

## Scope

This repository contains analysis code only. No raw data, prepared study data, sample identifiers, saved models, or study results are distributed. The code requires locally supplied input tables described in [data/README.md](data/README.md). Consequently, this repository alone does not reproduce the numerical results of the paper. No data-access policy or permission to redistribute study data is implied.

| Script | Analysis |
| --- | --- |
| `sliding_window_analysis.R` | Sample-level mixed models within windows; main Figure 2a-b and window-size sensitivity underlying S3. |
| `mediation.R` | Window-level statistical mediation underlying Figure 4a. |
| `lavaan_SEM_Fig3.R` | Separate lavaan SEMs for the Figure 3 groups and, if supplied, Figure S5 tertile groups. |
| `analysis_helpers.R` | Required shared validation, parallel execution and export functions. |
| `00_run_all.R` | Optional entry point; runs the three analyses and stops on an error. |

Data preparation from sequence counts, group assignment, interaction tests, random forests, the Figure S4 nested model, the Figure 4 SEM, and plotting are outside this selected release. The original local workflow's disabled window bootstrap/location-omission branches are not included. The publication release keeps the numerical methods of the included analyses while replacing machine-specific input discovery with documented relative input paths.

## Requirements and execution

The local study workflow used R 4.5.1. Exact package versions are not pinned in this release; every completed run writes `sessionInfo.txt`. Numerical equivalence across different package versions is not guaranteed.

Install dependencies yourself once:

```r
install.packages(c("lme4", "lmerTest", "mediation", "lavaan", "openxlsx"))
```

1. Supply the required UTF-8 CSV tables in `data/`; retain the original within-table row order.
2. Open a saved script in RStudio and click **Source**, or set the working directory to the repository root and use:

```r
source("sliding_window_analysis.R")
source("mediation.R")
source("lavaan_SEM_Fig3.R")
# Alternatively, after supplying every input:
source("00_run_all.R")
```

Do not supply Figure 2's full-rhizosphere sliding windows to the mediation script. Its input is a separately prepared paired-sample window table. Running the first script does not create that table or the SEM groups.

The mixed-model and SEM scripts default to two PSOCK workers for independent deterministic fits. Set `workers <- 1L` for serial execution. The mediation bootstrap stays serial to preserve its original seed and random-number sequence. No script installs packages, downloads data, or calls another analysis to silently generate missing inputs.

## Statistical specifications

### Sliding windows

The main model is:

```r
Multifunctionality ~ Multidiversity + (1 | Location) + (1 | Vegetation)
```

`Location` is sampling site; `Vegetation` is tree taxon, not individual tree identity. The original additional `Pathogen_adjusted` model adds log10 pathogen percentage as a fixed covariate. It is retained and labelled separately; it is **not** the nested random-effects sensitivity model used for Figure S4.

Windows of 50, 60, 70 and 80 samples advance by one sample. Samples are ordered by pathogen percentage and then SampleID to break ties. Fits use REML, `bobyqa`, Satterthwaite slope tests, and t-based 95% intervals. The standardized slope is `b * SD(Multidiversity) / SD(Multifunctionality)` using the full input soil compartment, not each window. The trend predictor is `mean(log10(pathogen percentage))`, not `log10(mean(pathogen percentage))`.

The window-60 main results are exported separately for RS and BS. Failed models remain flagged rows; singular but converged models are retained. OLS trends use valid window estimates. A run containing failed window fits is marked `NEEDS_REVIEW` and stops the combined entry point after exporting diagnostics.

### Mediation

For the published analysis, 304 paired samples give 245 windows of size 60 and step 1. Window slopes retain the standardization reference from all 306 rhizosphere samples. The input must already contain these prepared estimates.

The exposure X is the arithmetic-mean pathogen percentage, M the window-mean root-effect strength, and Y the standardized diversity slope. X is not log-transformed. X, M and Y are z-scored over complete windows before fitting `M ~ X` and `Y ~ X + M`. The script uses `mediation::mediate`, 5,000 ordinary bootstrap replicates, seed 2025, percentile 95% intervals, and a one-SD exposure contrast (0 to 1). Path P values come from OLS t-tests. The reported mediated proportion is ACME / total effect.

`mutiBC` (all fungal ASVs retained) is the primary analysis. If the separately prepared `mutiBC_NoPathogen` column is supplied, the original pathogen-excluded sensitivity case is also analyzed and explicitly labelled. Residual-correlation sensitivity is retained as a table when supported. Redundant alternative definitions of the mediated proportion and graphical sensitivity output are omitted from this selected code release.

Adjacent windows overlap. Ordinary OLS, SEM and window-level bootstrap inference do not correct that dependence. Mediation estimates describe observational associations and do not establish causal mediation.

### Group-specific SEM

The script contains the original six response equations and 26 directed paths. It uses `lavaan::sem`, ML, `meanstructure = TRUE`, and complete cases across all eight SEM variables. Each variable is z-scored within its **allocated** group with `na.rm = TRUE`, before incomplete observations are excluded. Do not prefilter the input or scale across groups.

Every reported path estimate, SE, P value and interval comes from `standardizedSolution(type = "std.all")`. The script does not output CommonScaleBeta or alternative coefficient scales. It reports lavaan response R-squared and global fit statistics; it does not perform a cross-group path-difference test or add random effects.

Group membership is supplied by the preceding grouping workflow. The primary groups follow the first adjacent significant-to-nonsignificant transition in the RS main size-60 windows; the high group is above the arithmetic-mean pathogen percentage of the first nonsignificant window, and the low group contains an equal number of the lowest-abundance RS samples. Do not replace this with a median split or a hard-coded rounded 6% cutoff. Tertile groups are optional. Group-specific sample counts before and after exclusions are saved.

## Outputs and data separation

Each run creates its own timestamped directory under `outputs/`, with CSV result tables, an audit workbook, `sessionInfo.txt`, and a status file. Outputs can contain private sample/window identifiers and standardized inputs. They are excluded by `.gitignore`, as are local data and common result formats. This does not prevent an explicit force-add or manual web upload: upload only the code and documentation supplied in this release, not your populated data/output folders.

No models were run while preparing this release. Files received static checks only; statistical runtime equivalence with the original workflow has not been tested. No demonstration data or fabricated results are included.
