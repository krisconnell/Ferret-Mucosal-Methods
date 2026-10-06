
Kristin Connelly, 3/2/2025 - Diversity and rarefaction analysis work 
**Work in progress**
Last Updated **10/1/2026** 

This folder contains scripts calculating B cell diversity across different tissues in a ferret model. It is in Haskell and in python and is designed to be run at the command line. The diversity analysis included is sourced from https://github.com/DrexelSystemsImmunologyLab/diversity. 

The scripts first take a search database as an input, collapse different clone reads with same CDRH3 sequences into smaller clone subsections, and outputs different fasta files that represent databases with x or more instances of each clonal lineage. This is helpful to weed out clones that are not sampled enough to be considered an accurate representation of the tissue environment. Note the input search database used for creating this script comes from Luke's new BCR pipeline using IgBLAST, which is present in the FerretBCR folder modified from https://github.com/LukeHebert/bcrseq_igseq. 

You will need to adjust your fasta headers first by adding a "sample/source" name to the beginning of  the entry: BAL|xxxxx, for example. This is completed in using the provided code and is described in the accompanying png. 

The first pipeline, run_clone_threshold_pipeline.py, is run assuming your directory (the animal) is separated into subdirectories for each tissue, which contains FASTA files for analysis, taken from the make_searchable.py function + rewritten naming system. It needs the following: 

"-i", "--input" - "Input FASTA file"
"-s", "--sample" - "Sample name"
"-t", "--thresholds" type=int, "Clone instance count thresholds"
	-I use [1, 5, 10, 15, 20, 30, 50]
"-o", "--output_prefix, "Output prefix for files"

These are already set in the scripts but adjust them according to your needs. 

The second part of the pipeline runs your resulting thresholded fasta files through the diversity tool found at https://github.com/DrexelSystemsImmunologyLab/diversity?tab=readme-ov-file. *check the Additional_Stack_info file for more information on setting this up.*

Needs: 

"-d", "--directory", "Directory containing thresholded FASTA files"
"-s", "--sample", "Sample name"
"-o", "--output_prefix", "Prefix for output diversity files"
"-n", type=int, default=40, "Number of subsamples (default: 40)"

Again, these are set in the scripts but need to be adjusted for your needs. 

Results in diversity calculations (Shannon, Simpson, hill) and rarefaction curves you can use to determine appropriate number of clones to use for tissue analysis. Can be used in downstream analysis as well for determining clonal overlap. 

<img width="2200" height="1540" alt="DiversityPipeline2" src="https://github.com/user-attachments/assets/d95eddcc-cbc2-495a-b64c-f9ddc8a54fa3" />
