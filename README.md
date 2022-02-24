# FANCONI

### STEPS in the analysis 

1. Analize each one of the fanconi samples separately
2. Annotate each of the with the label transfer function
3. Differential analysis corrected vs uncorrected per cell type
4. GSEA analysis in each one of the cell types


###INTEGRATION
1. Integrate each one of the samples with the healthy donor
2. Differential analysis healthy vs uncorrected and healthy vs corrected
3. GSEA analysis in the contrast healthy vs uncorrected

# SCRIPTS
 
01_Cellranger_counts.sh: script to create count matrix from fastq files

02_fanconi_analysis: all the steps of individual sample analysis (processing and annotation)
03_plots_fanconi_samples: all the code for the plots of the fanconi samples (figure 2)
04_GSEA: 
  
