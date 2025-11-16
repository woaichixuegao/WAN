# WAN: A Workflow for Invasive Species Habitat Preference Analysis

**WAN** is an R package designed to support a complete analytical workflow for invasive species habitat-preference assessment.
It provides tools for:

* Geospatial preprocessing
* Distance-to-habitat computation
* Frequency–distance decay analysis
* Trend modelling (LM and GAM)
* Statistical comparison among habitats
* Random-forest modelling combining habitat and climate predictors
* Reproducible case studies and pipeline utilities

The package integrates **sf**, **terra**, **ggplot2**, **mgcv**, and **randomForest** to provide an end-to-end framework for habitat-based ecological analysis.

---

## Key Features

### 1. Occurrence & Background Point Processing

* Import occurrence records from CSV
* Generate background points using:

  * Random sampling
  * Restricted buffers
  * Target-group background
  * Bias layers

### 2. Distance Computation

* Compute distance from each point to multiple habitat layers
* Supports `sf` and `terra` vector formats
* Multi-habitat distance extraction in one step

### 3. Distance–Frequency Analysis

* Build **long** and **wide** format frequency tables
* Automatic or manual distance-class breaks
* Exponential-decay fitting with confidence intervals

### 4. Trend Modelling

* Least-squares LM (linear or polynomial)
* GAM smoothing with `mgcv`
* 95% confidence intervals for fitted curves
* Multi-habitat trend comparison and summary tables

### 5. Statistical Testing

* Kolmogorov–Smirnov tests
* Permutation tests
* Decay-rate comparisons
* Bootstrap procedures with zero handling

### 6. Machine-Learning Integration

* Prepare presence–background data frames for ML
* Random-forest modelling
* Variable-importance ranking (habitat distances + climate variables)

### 7. One-Click Pipeline

`run_invasion_habitat_pipeline_generic()` performs the full workflow:

1. Read occurrence data
2. Generate background points
3. Compute distances
4. Build frequency tables
5. Produce analysis-ready objects

---

## Installation

### Development version (GitHub)

```r
install.packages("remotes")
remotes::install_github("woaichixuegao/WAN")
library(WAN)
```
---

## Basic Example

```r
library(WAN)

# Example input
occ_file <- "occurrence_points.csv"

habitats <- list(
  roads    = roads_sf,
  rivers   = rivers_sf,
  sand     = sand_sf,
  stations = stations_sf
)

res <- run_invasion_habitat_pipeline_generic(
  occ_file = occ_file,
  habitats = habitats,
  n_bg     = 5000
)

str(res$freq_long)
```

---

## Plotting

### Exponential Decay (All Habitats)

```r
plot_multi_habitat_decay(res$freq_long, ncol = 2)
```

### LM / GAM Trends

```r
plot_multi_habitat_trends(
  res$freq_long,
  methods = c("lm", "gam"),
  degree  = 1,
  k       = 5
)
```

---

## Random Forest Example

```r
ml_data <- prepare_ml_data(
  dist_occ = res$freq_long,
  dist_bg  = res$freq_long,
  env_occ  = climate_occ,
  env_bg   = climate_bg
)

rf_res <- fit_random_forest_importance(ml_data)

rf_res$importance
```

---

## Documentation

Full documentation is available in the package:

```r
help(package = "WAN")
```

A case-study vignette is included in the `vignettes/` folder of the repository.

---

## Project Structure

```
WAN/
│ DESCRIPTION
│ NAMESPACE
│ README.md
│
├─ R/                    # All R functions
├─ man/                  # roxygen2-generated documentation
├─ vignettes/            # Long-form documentation
├─ data-raw/             # Optional data scripts
└─ inst/                 # Additional resources
```

---

