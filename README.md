splitRtools: Preprocessing tools for SPLiT-seq data
================

[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/badge/devel%20version-0.0.1.1-blue.svg)](https://github.com/JamesOpz/splitRtools)
[![](https://img.shields.io/github/languages/code-size/JamesOpz/splitRtools.svg)](https://github.com/JamesOpz/splitRtools)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)

# Welcome to the splitRtools package!

## This package is under active development and all functionality is not yet validated!!

## The package may change significantly over development

## :arrow\_double\_down: Installation

The package can be installed from this github repository:

``` r
# Insall BiocManager if not present
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages 
BiocManager::install(c("zellkonverter", "ShortRead", "scater", "DropletUtils",
                        "ComplexHeatmap"))

# Install devtools for github installation if not present
require(devtools)

# Install package from github repo
devtools::install_github("https://github.com/JamesOpz/splitRtools")
```

## Overview

The splitRtools package is a collection of tools that are used to
process SPLiT-seq scRNA-seq data first published in [Rosenberg et.al,
2019](https://www.science.org/doi/10.1126/science.aam8999?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed).
</br> </br> The splitRtools package is designed to take as input data,
the various outputs from the [zUMIs
package](https://github.com/sdparekh/zUMIs)
([paper](https://academic.oup.com/gigascience/article/7/6/giy059/5005022?login=true))
for scRNA-seq barcode mapping and alignment. The zUMIs package is used
to take raw FASTQ output, assign and filter reads to barcodes, then map
the cDNA reads to a reference genome using STAR producing a CellxGene
matrix, as well as some reporting about the pipeline outputs. </br>
</br> A sample zUMIs pipeline with configuration to work with the
Rosenberg-2019 barcode setup is available
[here](https://github.com/JamesOpz/split_seq_zUMIs_pipeline).

## Running the splitRtools pipeline

### Data input directory structure

#### data\_folder

The splitRtools pipeline depends on the naming of the zUMIs pipeline
barcodes/read mapping output. All zUMIs outputs for each sublibrary must
be contained within a folder with the same name as the zUMI experiment
name. This is the name embedded into each zUMIs output file. The zUMIs
sublibrary output folder must also be named the same as this zUMIs
experiment name. The folders for each individual sublibrary must be
contained withing the `data_folder` and this folder’s absolute path must
be specified in the `run_split_pipe()` arguments. </br>

#### fastq\_path

The other input folder is the FASTQ folder containing the raw data used
as input for the zUMIs mapping pipeline. This allows zUMIs to calculate
the total reads from each sublibrary to calculate several metrics
relating to the experimental sequencing depth. The absolute path for
this folder is specified in the `fastq_path` arguments of the
`run_split_pipe()` function.</br>

#### File input structure

\|</br> \|–`data_folder`</br> \|          \|</br>
\|          \|-`sub_lib_1`</br>
\|          \|       \|-`sub_lib_1.BCstats.txt`</br>
\|          \|       \|-`zUMIs_output`</br> \|          \|</br>
\|          \|-`sub_lib_2`</br> \|          \|-`sub_lib_n`</br> \|</br>
\|-`fastq_path`</br>           \|</br>           \|`sub_lib_1`</br>
          \|       \|-`sub_lib_1_R1.fastq.gz`</br>
          \|       \|-`sub_lib_1_R2.fastq.gz`</br>           \|</br>
          \|-`sub_lib_2`</br>           \|-`sub_lib_n`</br> </br>

#### Barcode maps

The experiment barcoding layout must be provided as a csv file with two
columns - well position (numeric: 1-96) and barcode sequence in each
well. Currently `splitRtools` supports one barcoding layout for the RT
plate (args `rt_bc`) and another for the two subsequent ligation rounds
(args `lig_bc`). An example of the barcoding layout sheet (Rosenberg
2019 format) is located in this repository in `data/barcodes_v1.csv`.

#### Sample maps

Similar to the barcoding layout, the sample layout for the RT barcode
indexing needs to be provided as - well position and sample\_id. This
enables the labelling of each cell with its sample of origin and is
specified in arg `sample_map`. An example of the sample map layout sheet
is located in this repository in `data/cell_metadata.xlsx`.

### Executing the pipeline

The splitRtools pipeline is run through the `run_split_pipe()` function,
which acts as a wrapper function to execute the pipeline. A basic setup
for the pipeline is as follows: (for more information on pipeline
arguments use `?run_split_pipe`) </br>

``` r
# Load splitRtools
library(splitRtools)

# Run the splitRtool pipeline
# You must always point to two parent folders, 
# one containing zUMIs data folders and raw FASTQ data folders
run_split_pipe(mode = 'merge', # Merge sublibraries or process seperately
               n_sublibs = 2, # How many to sublibraries are present
               data_folder = "./../test_data_sp_5_miseq/", # Location of zUMIs data directory
               output_folder = "../test_data_sp_5_miseq_outputs/", # Output folder path
               filtering_mode = "knee", # Filter by knee (standard) or manual value (default 1000) transcripts
               fastq_path = "../fastq_single/", # Path to folder containing subibraru raw FASTQ
               rt_bc = "../test_data_sp_5_miseq/barcodes_v1.csv", # RT barcode map
               lig_bc = "../test_data_sp_5_miseq/barcodes_v1.csv", # Ligation barcode map
               sample_map = "../test_data_sp_5_miseq/cell_metadata.xlsx" # RT plate layout file
)
```

## Pipeline outputs

### Output directory structure

\|</br> \|–`output_folder`</br>           \|</br>
          \|-`sub_lib_1`</br>
          \|       \|-`unfiltered_sce_h5ad_objects`</br>
          \|       \|-`filtered_sce_h5ad_ojects`</br>
          \|       \|-`ggplot_outputs`</br>
          \|       \|-`report_data_outputs`</br>           \|</br>
          \|-`sub_lib_2`</br>           \|-`sub_lib_n`</br>
          \|-`merged_sublibrary_data`</br>

### Output data

The first stage of the pipeline labels converts the cell count matrix
into a `SingleCellExperiment` object and labels each cell with various
`ColData` with a series of well IDs based each stage of the barcoding
process and the correspondence between the RT wells ID and the
`sample_map` .xlsx file provided. This data is then stored as an `SCE`
or an `annData` object in `unfiltered/` output folder for each
sublibrary.</br> </br>

### Diagnostic plots

The splitRtools pipeline will generate a set of diagnostic plots in
order to evaluate the initial quality of the SPLiT-seq scRNA-seq data.
</br> </br> After labeling the data is filtered using either the
`DropletUtils` package spline-fitting functionality or a user specified
manual cutoff of transcripts. This produces the following waterfall plot
along with quantifiaction of the cell types recovered by sample: </br>
</br>
<img src="data/3_umi_waterfall.png" width="380"><img src="data/cell_abundance_barplot.png" width="200">
</br>  
</br> The barcoding cell data is then mapped to the respective plate
locations across the 3 barcoding rounds to provide a series of heatmaps
displaying cells recovered per well and median UMI per cell across all
wells: </br>
<img src="data/rt_barcoding_layout.png" width="400"><img src="data/rt_umi_layout.png" width="425">
