# Analysis Code for Wilcox et al. 2025

This repository contains the R code used for analyzing spatial, temporal, taxonomic, and phylogenetic bias patterns in biodiversity data sources across three California study regions (Santa Monica Mountains, Mount Tamalpais, and Marble/Salmon Mountains).

## Abstract
Digitized herbarium specimens and iNaturalist observations provide invaluable plant biodiversity data. Combining these two data sources could create a more holistic representation of local biodiversity; however, understanding biases inherent to each is critical to determine how to best combine and utilize these data. 
We examined how interpretation of taxonomic and phylogenetic diversity, naturalized species detection, and spatiotemporal coverage differ when using herbarium and iNaturalist data alone and together. We also examined how these patterns vary among areas with different degrees of collecting and community science efforts. 
Across areas, diversity was higher when data sources were combined, and complementary spatiotemporal coverage between data sources indicates that combining these data is useful; however, biases unique to each data source should be considered during analyses. Naturalized species detection, diversity patterns, and spatial biases varied by area indicating that local context impacts our current views of biodiversity and should shape future monitoring.
Our findings suggest that continued digitization and georeferencing of herbarium records will help provide critical information about biodiversity, but strategic collecting of both specimens and iNaturalist observations moving forward will ensure that we are capturing biodiversity change in real time, helping us track responses to environmental change.   

## Dataset
The analyses use a filtered and annotated GBIF dataset that is available at:

Hill, A. (2025). Filtered and Annotated GBIF dataset for Wilcox et al 2025 [Data set]. Zenodo. https://doi.org/10.5281/zenodo.15265921

## Analysis Scripts
The repository is organized into several R scripts that implement the analyses described in the manuscript:

### Data Preparation
- `00-data-prep.R`: Processes GBIF plant occurrence data for the three study regions, joining with boundary data and verifying taxa against a curated species list.
- `00-trail-prep.R`: Extracts and processes road and highway data from OpenStreetMap, applying a 3km buffer around each study area.

### Spatial Analysis
- `01-NNI-analysis.R`: Implements Nearest Neighbor Index analysis to quantify spatial clustering patterns between iNaturalist observations and herbarium specimens.
- `01-study-area-map.R`: Generates visualizations of the three study areas showing their locations, ecoregions, and elevation ranges.
- `02-nearest-road.R`: Calculates distance from each occurrence point to the nearest road or trail using OpenStreetMap data.
- `02.5-road-length-and-density.R`: Calculates total road length and road density (km/km²) for each study area.
- `02.5-road-regression-analysis.R`: Performs Gamma regression models to analyze the relationship between occurrence locations and distance to roads.
- `03-spatial-bias-figure.R`: Creates the multi-panel figure showing spatial bias patterns, including road networks, occurrence density, and distance distributions.

### Temporal Analysis
- `01-Temporal_coverage_bias.R`: Analyzes and visualizes species accumulation through time, native vs. non-native species accumulation patterns, and changes in species richness, record numbers, and collector/observer numbers over time.

### Taxonomic and Phylogenetic Analysis
- `01-Taxonomic_coverage_bias.R`: Examines patterns of taxonomic bias between data sources, analyzing family-level representation and specialization.
- `01-Phylogenetic_coverage_bias.R`: Implements phylogenetic analyses to assess evolutionary bias in biodiversity data.

### Miscellaneous
- `03-full-dataset-n-citation.R`: Prepares dataset for uploading to Zenodo.

## Usage Notes
 The required R packages include:
- Spatial analysis: `sf`, `stars`, `spatialEco`
- Data manipulation: `tidyverse`, `duckdb`
- Visualization: `ggplot2`, `patchwork`, `cowplot`
- Stats: `mgcv` (for GAM models)

