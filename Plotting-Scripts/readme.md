# Analyzing thresholded animal and tissue-specific data
First, this set of scripts is meant to analyze the diversity outputs from run_diversity_batch.py using the analyze_diversity_outputs.py script, which graphs diversities across selected thresholds for your total animal and for specific tissues in your animal as a function of retained mean and median diversity over increasingly stringent thresholding. These graphs also include separate tissue-specific diversity graphs and also report a reccomended threshold for further analysis which considers global and tissue-specicic retained diversity in calculations (i.e. the elbow in the elbow plot, typically). 

After getting a thresholded value, apply_threshold_and_merge_global.py applies the threshold to tissue datasets and global datasets, and it also creates IgG, IgA, and IgM-specific datasets globally and at the tissue-specific level. 

# Scripts for repertoire analysis: 

- **chord_diagram_tissues.py** : use Python holoviews package to create a chord diagram representing lineage sharing across tissues in an animal, as chords connecting two separate tissue sources. Note that this does not visualize 3+ tissue connections.
- **repertoire_overlap_metrics.py** : Measure similarities across tissues via cosine similarity, jaccard similarity, morisita-horn similarity, lineage overlap, and lineage overlap as percent abundance of a repertoire (preserving directionality).
- **CDRH3LengthPlot.py** : creates a PNG file of graphed CDRH3 lengths across tissues with mean values indicated.
- **SHMPercentbytissue.py** : creates a PNG file of graphed SHM across tissues with mean values indicated. SHM calculated as difference in v_identity percent.
- **isotype_overlap.py** : creates PNG files of IgG, IgA, and IgM lineage overlap across all tissues. Also has an option to define tissue "groups" like LRT, URT, systemic, et cetera, editable within the code.
- **tissue_hydropathy_plot.py** : creates a PNG file gauguing repertoire hydropathy across tissues.
- **Get logos plotting script from Luke and add it here!**
  
<img width="3000" height="2100" alt="PlottingScripts" src="https://github.com/user-attachments/assets/189eb700-1077-46f2-92a9-7433a63ac3fa" />

