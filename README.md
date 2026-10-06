**DRAFT** 

# Ferret-Mucosal-Methods
Repository containing scripts for reading BCR fastq files, identifying genes via IgBLAST, clustering lineages at the global level, splitting repertoires into tisssue-specific responses, and completing downstream diversity analysis using code from https://github.com/DrexelSystemsImmunologyLab/diversity.

The original BCRseq pipeline is created by https://github.com/LukeHebert with minimal edits for use in this pipeline. Some notable changes: 
- clustering is baseed on 90% CDRH3 identity here using cluster.py compared to Luke's gupta_cluster.py, which groups sequences by identical V gene, J gene, and CDRH3 amino acid (cdr3_aa) length.
- the animal repertoire and tissue-specific repertoire undergoes thresholding to account for different sampling depths across tissues.
- A series of scripts are included for plotting data after analysis. 
