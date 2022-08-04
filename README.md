# PD-cell-type-MR
Codes used in PD MR articles

#Exposure data
Expression quantitative trait loci (eQTL) data of dopaminergic neurons from a recent study by Jerber et al.1 were used in this study. 
In brief, Jerber et al. differentiated 215 human-induced pluripotent stem cell (iPSC) lines (from the HipSci consortium2) into several neuronal cell types, 
including dopaminergic neurons. Single-cell RNA sequencing was performed to obtain gene expression levels in each cell type. 
To investigate the associations between common genetic variations and gene expression in each differentiated cell type, 
expression quantitative trait loci analyses were conducted by integrating genotype data (from the HipSci consortium) and the single-cell RNA sequencing results.
Cells were harvested at several differentiation time points (11, 30, and 52 days) for single-cell RNA sequencing, 
and only eQTL data from day 52 (differentiation of iPSC toward dopaminergic neurons for 52 days) were used in this study. 
To investigate the effect of oxidative stress on gene expression, 
Jerber et al. also conducted an eQTL analysis by treating the differentiated dopaminergic neurons with rotenone. 
More details about the eQTL of dopaminergic neurons can be found in the study of Jerber et al.

#Outcome data
The largest PD GWAS reported by Nall et al.3 was used in this study. Briefly, Nall et al. 
performed a large-scale GWAS meta-analysis on PD by combining data from 17 cohorts, 
including 37,688 PD cases, 18,618 UK Biobank proxy-cases and about 1.4 million controls. 
As GWAS results of several datasets (such as 23andme) were not publicly available, 
genome-wide associations of 33,647 cases and 449,056 controls (referred as main GWAS) were used in this study. 
We noticed that some SNPs were not available in the main GWAS. 
However, association results of these SNPs were available in a relatively smaller PD GWAS (26,421 cases and 442,271 controls, also from by Nall et al ). 
To include more variants, for SNPs that were not available in main GWAS, association results from 26,421 PD cases and 442,271 controls were used. 
Detailed information about the PD GWAS can be found in the original publication3 and sample size information for each IVs can be found in Supp_tables.xlsx, 
“IVs_sample_size” sheet.

#Instrument selection
Due to the small sample size of the eQTL dataset used (N=215 iPSC lines), 
to include a proper number of genetic instrumental variables (IVs) for MR analysis and to ensure that the basic assumptions of MR are satisfied, 
i.e., the instrument variable has a strong impact on exposure (SNPs are strongly correlated with gene expression)4, 
we used a relatively relaxed eQTL P threshold, as described in previous studies 5-7. 
These studies showed that filtering eQTL data with 1 × 10−5 guarantees a small weak instrument bias. 
For the above reasons, we chose a less stringent threshold (cis-eQTL P Val < 1×10-5) as the cis-eQTL threshold in this study. 
To ensure the reliability of MR results, we adopted a strict significance threshold 
(Bonferroni correction threshold) to control the potential false positive results of MR analysis.
For each exposure dataset, the SNPs that were associated with gene expression (cis-eQTL P < 1×10-5) 
in dopaminergic neurons were selected as genetic instrumental variables (IVs). 
To obtain independent IVs, we firstly conducted LD clumping using the PLINK clumping method. 
Genotype data of Europeans from the 1000 Genomes were used as a reference panel. 
In a 10 Mb window, if LD values (r2) of two or more SNPs were smaller than 0.001, these SNPs were considered independent IVs.

#Sensitivity analyses
For exposures with multiple IVs, we also tested the heterogeneity across variant-level MR estimates, 
using the “mr_heterogeneity()” function in TwoSampleMR package (Cochrane Q method), 
and removed exposures with significant heterogeneity (Q_pval < 0.05). 
In addition, to test whether there is horizontal pleiotropy among multiple IVs (more than 2 IVs), 
a pleiotropy test was also performed by using MR Egger analysis. 
We tested the difference between the intercept term and 0, and the exposure with statistical difference (P-val > 0.05) was removed. 
Finally, if there are exposures with more than 2 IVs, a Leave-one-out sensitivity test is also performed, 
and the MR results of the remaining IVs are calculated by removing the IVs one by one to ensure the robustness of the MR results. 
Please refer to Supp_tables.xlsx for the MR results and sensitivity analysis results.

#Mendelian randomization analysis
To ensure the effectiveness of an SNP on exposure and the effect of that SNP on outcome each corresponds to the same allele, 
exposure and outcome data were harmonized using the “harmonise_data()" function from the TwoSampleMR package8. 
When only one SNP was available as a genetic instrument, we used the Wald ratio test for MR analysis, 
i.e. the ratio between the instrument-outcome and instrument-exposure association estimates. 
If there were two or more IVs, the inverse variance weighted (IVW) approach was used. 
The weighted sum of the contribution of each independent IV constitutes the final results of the IVW approach. 
The MR results were corrected by Bonferroni multiple correction approach. For the dopaminergic neuron-specific eQTL data, 
a total of 1,028 genes were included in the final MR analysis, so the Bonferroni correction threshold was 4.86 x 10-5. 
For the rotenone-treated dopaminergic neurons eQTL data, a total of 640 genes were included in the final MR analysis, 
so the Bonferroni correction threshold was 7.81 x 10-5. Genes passed the significance threshold were considered to be significant. 
