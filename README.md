[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub](https://img.shields.io/github/license/Monashbioinformaticsplatform/LFQ-Analyst?color=brightgreen)
![R](https://img.shields.io/badge/R-%3E4.2-brightgreen)

# LFQ-Analyst
A tool for analysing label-free quantitative proteomics dataset https://bioinformatics.erc.monash.edu/apps/LFQ-Analyst/

![LFQ-analyst_pipeline](./www/static/img/LFQ_analyst.svg)






## Motivation

- Automate downstream statistical analysis of Label free quantitative proteomics data (generated by MaxQuant)


### Input

- MaxQuant **proteinGroups.txt** file
- An experiment design table (tab separated file) containing three columns ("label", "condition", "replicate")

#### Data pre-filtering criteria

- Remove potential contaminants
- Remove reverse sequences
- Remove proteins identified only by sites
- Remove proteins identified/quantified by a single Razor or unique peptide
- Remove observation with high proportion of missing values (intensity values must be present
at least 2 out of three replicates)

#### Advanced parameters to choose

- Differencial expression cutoff
  - Adjusted p-value cutoff (FDR cutoff on quantitation)
  - Log2 fold change cutoff
- Option to choose paired test for matched pair data
- Types of imputation
  - A number of missing value imputation options including knn, Minpob etc.
- Type of FDR correction
  -   Benjamin Hochberg (BH) method
  -   t-statistics correction: Implemented in
    [fdrtools](http://strimmerlab.org/software/fdrtool/)
- Option to include proteins identified/quantified with a single unique peptide.
- Select how many clusters of differentially expressed proteins needed for the heatmap (default is 6)



### Outputs

#### Result table

-   **LFQ Results Table:** Includes names (Gene names), Protein Ids, Log
    fold changes/ ratios (each pairwise comparisons), Adjusted
    *p-values* (applying FDR corrections), *p-values*, Boolean values
    for significance, average protein intensity (log transformed) in
    each sample.

#### Result Plots
  1. Interactive volcano plot for each pairwise comparison.
  2. Heatmap of differencially expressed proteins
  3. Protein intensity plots for a single or group of selected proteins from table. 

#### QC Plots
  1. PCA plot (Could move to QC section)
  2. Sample Correlation (pearson correlation)
  3. Sample Coefficient of variations (CVs)
  4. Number of proteins per sample
  5. Sample coverage (overlap of identified proteins across every sample)
  6. Missing value heatmap
  7. Imputation effect on sample distribution

### Download options

**Download tables** (csv format)

1.  Results: Same as *LFQ Results Table*
2.  Unimputed data matrix: Original protein intensities before
    imputation in each sample.
3.  Imputed data matrix: Protein intensities after performing selected
    imputation method
4.  Full results: Combined table of all above data outputs i.e. with and
    without imputation information, along with fold change and p-values.

**Download Report** 
- A summary report for each analysis that
    includes method, summary statistics and plots.


### Local installation

The current version of LFQ-Analyst is hosted on `R - 4.2.1`. The detailed dependency information can be found in the `dependencies.txt` file.

Once installed all the dependencies following steps to run the server locally.

- Using git and Rstudio
```
## Clone the repository
git clone https://github.com/MonashBioinformaticsPlatform/LFQ-Analyst.git

## Move to the folder
cd LFQ-Analyst

## Inside R console or R studio
> library("shiny")

> runApp()

```

- Using Docker

Install & start Docker demon on your PC

```
## Option one：
## Pull LFQ-Analyst image from Docker Hub (From terminal)
> docker pull haileyzhang/lfq-analyst:tagname

## Option two：
## Clone the repository
git clone https://github.com/MonashBioinformaticsPlatform/LFQ-Analyst.git

## Move to the folder
cd LFQ-Analyst

## Build LFQ-Analyst (Any name after -t)
> docker build -f Dockerfile -t LFQ-Analyst .

## Run LFQ-Analyst (From terminal)

> docker run -p 3838:3838 LFQ-Analyst

## Open local interface

https://localhost:3838/LFQ-Analyst


```
