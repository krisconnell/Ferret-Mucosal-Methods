Kristin Connelly, 3/2/2025 - Diversity and rarefaction analysis work 
**Work in progress**

This folder contains scripts and example data for calculating B cell diversity across different tissue types in a ferret model. It is in Haskell and python and is designed to be run at the command line. 

The scripts first take a search database as an input, collapse different clone reads with same CDRH3 sequences into smaller clone subsections, and outputs different fasta files that represent databases with x or more instances of each clonal lineage. This is helpful to weed out clones that are not sampled enough to be considered an accurate representation of the tissue environment. Note the input search database used for creating this script comes from Luke's new BCR pipeline using IgBLAST: Z:\kmc5996\Ferret_Mucosal_Immunity\bcr_analysis_scripts_2024-09-25A

you will need to adjust your fasta headers first by adding a "sample/source" name to the beginning of  the entry: BAL|xxxxx, for example. 

After rewriting, start with the next part of the pipeline:

The first pipeline, run_clone_threshold_pipeline.py, is run assuming your directory (the animal) is separated into subdirectories for each tissue, which contains FASTA files for analysis, taken from the make_searchable.py function + rewritten naming system. It needs the following: 

"-i", "--input" - "Input FASTA file"
"-s", "--sample" - "Sample name"
"-t", "--thresholds" type=int, "Clone instance count thresholds"
	-I use [1, 5, 10, 15, 20, 30, 50]
"-o", "--output_prefix, "Output prefix for files"

The second part of the pipeline runs your resulting thresholded fasta files through the diversity tool found at https://github.com/DrexelSystemsImmunologyLab/diversity?tab=readme-ov-file. Needs: 

"-d", "--directory", "Directory containing thresholded FASTA files"
"-s", "--sample", "Sample name"
"-o", "--output_prefix", "Prefix for output diversity files"
"-n", type=int, default=40, "Number of subsamples (default: 40)"

Results in diversity calculations (Shannon, Simpson, hill) and rarefaction curves you can use to determine appropriate number of clones to use for tissue analysis. Can be used in downstream analysis as well for determining clonal overlap. 

NEXT STEPS: make a script that takes files and extracts thresholded subclass and expected richness values as well as diversity calcs. Keep the diversity number per threshold and record it. Then normalize the subclass and expected richness values and graph each thresholded limit all at the same time to gauge how quickly sampling gets to an appropriate level. Also report the threshold where we sufficiently reach 90% of total diversity with a low degree of subsampling (~10% of total sample? Determine this later.)


Make sure the following are in your path (at least where it is stored in my POD space): 

export PATH="$HOME/bcelldiversity_env/bin:$PATH" 
export PATH="/stor/work/Georgiou/kmc5996/Scripts/Bcelldiversity/diversity/.stack-work/install/x86_64-linux-tinfo6/2e42b64e23ad59927bc3d9de54d617fad19ea72f747c90575494cf6b45653aab/8.4.3/bin:$PATH"
export PATH="/usr/bin:$PATH"

Run this script using tmux! 


