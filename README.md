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

04_GSEA: script for gene set enrichment analysis in each one of the fanconi samples and code for plotting the ridgeplots with our own function our_ridgeplot.R also attached in the folder

05_integration.sh: script to run integration in the server

05_integration_fanconi_healthy: integration of one fanconi sample with the healthy

05_integration_SP: integration of peripheral blood samples

06_integration_DEA: integrated sample analysis and differential expression analysis of these samples

AUCell: explanation of aucell analysis with our signature and some plots

binomial: script to calculate the binomial of each sample

GSEO_microarrays: small analysis to compare our data and some public microarray data

healthy_MO268: script of the processing of the healthy sample then used in all the analysis

heatmap_rowannotation: script to plot the heatmap with row annotation

integration_fanconis: script to integrate all the fanconi samples

MYC_p53 : script with GSEA analysis of these two pathways myc and p53. It was done to compare with the public paper

our.likert: function to plot number of up and down regulated genes. Used in script 03

our_ridgeplot: function to plot our ridgeplot with GSEA results, used in script 04

proportions: script to calculate proportion test

wilcoxon: script to do wilcoxon and annova analysis to compare fanca expression levels





  
