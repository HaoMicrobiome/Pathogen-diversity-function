# Input specifications (no data included)

Supply your own UTF-8 CSV files locally. Empty fields and `NA` denote missing values. Do not upload study tables to this code-only repository. Column names are case-sensitive. This file documents schemas, not a downloadable dataset.

## 1. `sliding_samples_RS.csv` and `sliding_samples_BS.csv`

One row per retained sample, before any windowing. RS is rhizosphere and BS is bulk soil. The manuscript inputs have 306 RS and 304 BS samples.

| Column | Type / meaning |
| --- | --- |
| SampleID | Unique nonempty sample identifier; used to break abundance ties. |
| Location | Sampling-site grouping identifier. |
| Vegetation | Tree-taxon grouping identifier (not an individual tree). |
| Multidiversity | Mean of bacterial and fungal richness, each independently min-max normalized within the full soil compartment. |
| Multifunctionality | Mean of CB, BG, NAG, LAP and AP activities, each independently min-max normalized within the full soil compartment. |
| PlantPathogenRelativeAbundance | Percentage of all fungal reads assigned to the existing putative plant-pathogen ASV set; strictly positive, at most 100. |

Indices must be prepared before windowing. Do not normalize separately within windows. Missing model variables are not silently dropped. Raw ASV counts, guild assignment and enzyme measurements are not read by this publication script; data preparation is outside this release.

## 2. `mediation_windows.csv`

ASVs are amplicon sequence variants. The exposure is pathogen abundance, the mediator is root-effect strength, and the response is the diversity effect.

One row per 60-sample window rebuilt from the **304 matched RS/BS pairs**, in ascending WindowID order. With step 1, there are 245 windows, not the 247 windows obtained from all 306 RS samples in Figure 2. Retain source row order for exact bootstrap reproducibility.

| Column | Type / meaning |
| --- | --- |
| WindowID | Unique sequential numeric window identifier. |
| WindowSize | Always 60 for this analysis. |
| PlantPathogenRelativeAbundance | Arithmetic mean of pathogen percentages in the paired-sample window; not a log mean. |
| root_effect_strength | Window-mean root-effect strength: first calculate bacterial and fungal Bray-Curtis dissimilarities between each paired rhizosphere and bulk soil sample using ASV relative abundances; average the two dissimilarities with equal weights, then average across pairs in the window. All fungal ASVs are retained. |
| DiversityEffect | Window main-model slope, multiplied by SD(Multidiversity)/SD(Multifunctionality) calculated from the full 306-sample RS reference. |
| root_effect_strength_without_pathogens | Optional window-mean root-effect strength for sensitivity analysis. Use the same calculation as above, but remove plant-pathogenic fungal ASVs from both soil compartments and renormalize the remaining fungal relative abundances before computing fungal Bray-Curtis dissimilarity. The bacterial component is unchanged. |

Root-effect strength is between 0 and 1. Keep exposure and mediator on their pre-z-score scales; the script standardizes across complete windows. This input cannot be generated merely by renaming columns in the Figure 2 sliding-window export. Use only valid matched-window estimates; do not fabricate a window mediator from sample-level or unmatched tables.

## 3. `sem_groups.csv`

One row per **allocated** sample per scheme. Required columns:

| Column | Type / meaning |
| --- | --- |
| SampleID | Unique within each scheme; membership must not overlap Low and High in the same scheme. Reuse across schemes is expected. |
| Scheme | `Threshold_main` (required); `TopBottom102_sensitivity` (optional). |
| PathogenGroup | `Low` or `High`; both must be supplied for every scheme. |
| MAT | Mean annual temperature (degrees Celsius). |
| MAP | Mean annual precipitation (mm). |
| pH | Soil pH. |
| TC | Soil total carbon in the original consistent measurement unit. |
| EctomycorrhizaeRelativeAbundance | EcM fungal percentage of total fungal reads. |
| PlantPathogenRelativeAbundance | Putative plant-pathogenic fungal percentage of total fungal reads. |
| Multidiversity | Existing RS multidiversity index, before group Z standardization. |
| Multifunctionality | Existing RS nutrient-acquisition function index, before group Z standardization. |

Retain allocated samples with missing SEM covariates as NA: the original analysis computes each variable's mean/SD using its available values within the allocated group, then fits common complete cases. The manuscript's initial groups contain 86/86 samples for the primary scheme and 102/102 for the tertile scheme. Complete-case counts are 86/85 and 102/101, respectively. The script reports the actual supplied counts rather than forcing these values.

This release consumes group membership; it does not recalculate the threshold or assign groups. The primary threshold comes from the first adjacent P < 0.05 to P >= 0.05 transition of the RS main size-60 sliding-window slope, taking the latter window's arithmetic-mean pathogen percentage. Preserve the original allocation and tie-breaking decisions.
