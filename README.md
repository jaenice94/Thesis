# Thesis

This repository contains scripts with which the data analysis for "Decoding Gene Regulatory Elements in Parkinson's Disease-Related Cell Types: A Genomic Approach" were performed. H3K27ac, H3K4me3, FOXA1-CUT&Tag data, ATAC-seq data, mouse snRNA-seq/RNA-seq data, promoter-capture HiC data used for this analysis were generated in the Nott lab and are part of UPDATE ONCE PUBLISHED. 

Preprocessing was performed first using nf-core/cutandrun, nf-core/atacseq, nf-core/rnaseq pipelines, followed by additional processing using scripts provided under preprocessing/. Peak calling was optimised by source of input (e.g. sorted nuclei or cells), assay (atac-seq/CUT&Tag) or target antibody (H3K27ac/transcription factor). 

Downstream analysis included
- visualisation of genome-browser tracks (genome_tracks/)
- differential analysis (differential_analysis/)
- cell type and pathway enrichment analysis (pathway_and_enrichment/)
- transcription factor motif and activity enrichment analysis (tf_motif_analysis/)
- analysis of risk variant enrichment in cell type gene regulatory regions (ldsc/), including subanalysis to species conservation (conservation/)
- prioritisation and reannotation of PD risk genes (hmagma/)


Analysis tools used can be found in the scripts, details of the methods with versions and references will be made available once published UPDATE ONCE PUBLISHED. 
