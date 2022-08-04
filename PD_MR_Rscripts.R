
library("TwoSampleMR")
test_eqtl_exp_dat <-read_exposure_data("final_filter_1e-05_D52.DA.ROT_treated.qtl_results.txt",
                                       snp_col="chrsnp.id",beta_col="beta",se_col="beta_se",
                                       effect_allele_col="assessed_allele",
                                       other_allele_col="other_allele",
                                       pval_col="p_value",gene_col="feature_id",
                                       phenotype_col = "feature_id",
                                       sep = "\t", clump=TRUE)
write.table(test_eqtl_exp_dat,file="filter_1e-05_D52.DA.ROT_treated_MR_exp_dat.txt", row.names=F, quote=F,
            col.names=T, sep="\t")
test_eqtl_exp_dat <- read.table("ROT_treated/filter_1e-05_D52.DA.ROT_treated_MR_exp_dat.txt",header = T)
head(test_eqtl_exp_dat)

test_eqtl_exp_dat$id.exposure <- test_eqtl_exp_dat$SNP

outcome<-read.table("ROT_treated/overlap_PD_GWAS_D52.DA.ROT_treated",header = T)
head(outcome)
outcome$samplesize <- outcome$N_cases + outcome$N_controls
write.table(outcome,file="final_overlap_PD_GWAS_D52.DA.ROT_treated", row.names=F, quote=F,
            col.names=T, sep="\t")

outcome_dat <- read_outcome_data(
  filename = "ROT_treated/final_overlap_PD_GWAS_D52.DA.ROT_treated",sep = "\t",snp_col = "SNP",
  beta_col = "Beta",se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2",
  pval_col = "P",samplesize_col = "samplesize"
)

outcome_dat <- read_outcome_data(
  filename = "final_D52.DA.ROT_treated_nallsEtAl2019_TwosampleMR.txt",sep = "\t",snp_col = "SNP",
  beta_col = "Beta",se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2",
  pval_col = "P",samplesize_col = "samplesize"
)


dat <- harmonise_data(test_eqtl_exp_dat, outcome_dat)


##heterogeneity test
het <- mr_heterogeneity(dat)
het

#Q_pval Bonferroni threshold
0.05/7
#0.007142857

het_fialed <- het %>% dplyr::filter(Q_pval < 0.007142857)

## If there is strong heterogeneity among IVs (Q_pval is less than Bonferroni's correction threshold), then we need to eliminate these IVs

het_fialed
## There are no heterogeneous results



##Pleiotropy test
pleio <- mr_pleiotropy_test(dat)
pleio

# If there is no statistical difference between Egger_Intercept and 0 of MR-Egger (PVAL > 0.05), we can assume that there is no horizontal pleiotropy.

###Leave-one-out sensitivity test
single <- mr_leaveoneout(dat)
p <- mr_leaveoneout_plot(single)
p[["APIP.01dGB2"]]

p1 <- mr_scatter_plot(res, dat)
p1[["APIP.01dGB2"]]

#Single SNP analysis
res_single <- mr_singlesnp(dat)
res_single

p2 <- mr_forest_plot(res_single)
p2[["APIP.01dGB2"]]

p4 <- mr_funnel_plot(res_single)
p4[["APIP.01dGB2"]]



###MR analysis was continued

res <- mr(dat)
res <- subset_on_method(res)

0.05/640
#7.8125e-05
write.table(res,file="D52.DA_ROT_treated_MR_res.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

write.table(het,file="D52.DA_ROT_treated_MR_het_test_res.txt", row.names=F, quote=F,
            col.names=T, sep="\t")
            
     
##plot     
bleeding <- read.table("anno_D52.DA_ROT_treated_MR_res.txt", head=TRUE, sep = "\t")

library(qqman)
require("RColorBrewer")
display.brewer.all(type = "qual")

manhattan(bleeding,chr="feature_chromosome",bp="feature_start",p="pval",snp="exposure",
          col =c("#ec4035", "#1799d1"),
          main = "D52.DA_ROT_treated",genomewideline =-log10(7.8125e-05),
          suggestiveline = F,annotatePval = 7.8125e-05)





##########D52.DA_untreated##############


test_eqtl_exp_dat <-read_exposure_data("final_filter_1e-05_D52.DA.untreated.qtl_results.txt",
                                       snp_col="chrsnp.id",beta_col="beta",se_col="beta_se",
                                       effect_allele_col="assessed_allele",
                                       other_allele_col="other_allele",
                                       pval_col="p_value",gene_col="feature_id",
                                       phenotype_col = "feature_id",
                                       sep = "\t", clump=TRUE)
write.table(test_eqtl_exp_dat,file="filter_1e-05_D52.DA.untreated_MR_exp_dat.txt", row.names=F, quote=F,
            col.names=T, sep="\t")
test_eqtl_exp_dat <- read.table("D52.DA_untreated/filter_1e-05_D52.DA.untreated_MR_exp_dat.txt",header = T)
head(test_eqtl_exp_dat)

test_eqtl_exp_dat$id.exposure <- test_eqtl_exp_dat$SNP

outcome<-read.table("D52.DA_untreated/overlap_GWAS_D52.DA.untreated_for_MR.txt",header = T)
head(outcome)
outcome$samplesize <- outcome$N_cases + outcome$N_controls
write.table(outcome,file="D52.DA_untreated/final_overlap_PD_GWAS_D52.DA.untreated", row.names=F, quote=F,
            col.names=T, sep="\t")

outcome_dat <- read_outcome_data(
  filename = "D52.DA_untreated/final_overlap_PD_GWAS_D52.DA.untreated",sep = "\t",snp_col = "SNP",
  beta_col = "Beta",se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2",
  pval_col = "P",samplesize_col = "samplesize"
)



dat <- harmonise_data(test_eqtl_exp_dat, outcome_dat)


##Heterogeneity test
het <- mr_heterogeneity(dat)
het

write.table(het,file="final_overlap_PD_GWAS_D52.DA.untreated_het_result.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

#Q_pval Bonferroni threshold
0.05/6
#0.008333333

het_fialed <- het %>% dplyr::filter(Q_pval < 0.008333333)



het_fialed
## There are no results have heterogeneous 



##Pleiotropy test
pleio <- mr_pleiotropy_test(dat)
pleio


###Leave-one-out sensitivity test
single <- mr_leaveoneout(dat)
p <- mr_leaveoneout_plot(single)
p[["C7orf50.98iEy4"]]

p1 <- mr_scatter_plot(res, dat)
p1[["C7orf50.98iEy4"]]

#Single SNP analysis
res_single <- mr_singlesnp(dat)
res_single

p2 <- mr_forest_plot(res_single)
p2[["NDUFAF2.98iEy4"]]

p4 <- mr_funnel_plot(res_single)
p4[["C7orf50.98iEy4"]]



###MR analysis was continued

res <- mr(dat)
res <- subset_on_method(res)

0.05/1028
#4.863813e-05
write.table(res,file="D52.DA_untreated_MR_res.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

write.table(het,file="D52.DA_untreated_MR_het_test_res.txt", row.names=F, quote=F,
            col.names=T, sep="\t")
            
##plot
bleeding <- read.table("anno_D52.DA_untreated_MR_res.txt", head=TRUE, sep = "\t")

library(qqman)
require("RColorBrewer")
display.brewer.all(type = "qual")

manhattan(bleeding,chr="feature_chromosome",bp="feature_start",p="pval",snp="exposure",
          col =c("#ec4035", "#1799d1"),
          main = "D52.DA_untreated",genomewideline =-log10(4.863813e-05),
          suggestiveline = F,annotatePval = 4.863813e-05)











###########microglial###########

### SVZ Microglial #####

getwd()


library("TwoSampleMR")
test_eqtl_exp_dat <-read_exposure_data("anno_filter_1e-05_SVZ_eur_expression_peer5.cis_qtl_nominal.txt",
                                       snp_col="variant_id",beta_col="slope",se_col="slope_se",
                                       effect_allele_col="effect_allele",
                                       other_allele_col="other_allele",
                                       pval_col="pval_nominal",gene_col="SYMBOL_id",
                                       phenotype_col = "SYMBOL_id",
                                       sep = "\t", clump=TRUE)

write.table(test_eqtl_exp_dat,file="filter_1e-05_SVZ_eur_expression_exposure_data.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

test_eqtl_exp_dat <- read.table("microglial/filter_1e-05_SVZ_eur_expression_exposure_data.txt",header = T)
head(test_eqtl_exp_dat)

test_eqtl_exp_dat$id.exposure <- test_eqtl_exp_dat$SNP

outcome<-read.table("microglial/overlap_eQTL_nallsEtAl2019.txt",header = T)
head(outcome)
outcome$samplesize <- outcome$N_cases + outcome$N_controls
write.table(outcome,file="microglial/final_overlap_eQTL_nallsEtAl2019.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

outcome_dat <- read_outcome_data(
  filename = "microglial/final_overlap_eQTL_nallsEtAl2019.txt",sep = "\t",snp_col = "SNP",
  beta_col = "Beta",se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2",
  pval_col = "P",samplesize_col = "samplesize"
)



dat <- harmonise_data(test_eqtl_exp_dat, outcome_dat)


##Heterogeneity test
het <- mr_heterogeneity(dat)
het

## There are no results have heterogeneous 


##Pleiotropy test
pleio <- mr_pleiotropy_test(dat)
pleio



###keep doing MR analysis

res <- mr(dat)
res <- subset_on_method(res)

0.05/237
#0.0002109705
write.table(res,file="Microglial_PD_MR_res.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

##plot

bleeding <- read.table("microglial/anno_Microglial_PD_MR_res.txt", head=TRUE, sep = "\t")
head(bleeding)
library(qqman)
require("RColorBrewer")
display.brewer.all(type = "qual")

manhattan(bleeding,chr="CHROM",bp="POS",p="pval",snp="exposure",
          col =c("#ec4035", "#1799d1"),
          main = "Microglial_PD_MR_res",genomewideline =-log10(0.0002109705),
          suggestiveline = F,annotatePval = 0.0002109705)



### THA Microglial #####

test_eqtl_exp_dat <-read_exposure_data("anno_filter_1e-05_THA_cis_eQTL.txt",
                                       snp_col="variant_id",beta_col="slope",se_col="slope_se",
                                       effect_allele_col="effect_allele",
                                       other_allele_col="other_allele",
                                       pval_col="pval_nominal",gene_col="SYMBOL_id",
                                       phenotype_col = "SYMBOL_id",
                                       sep = "\t", clump=TRUE)

write.table(test_eqtl_exp_dat,file="filter_1e-05_THA_eur_expression_exposure_data.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

test_eqtl_exp_dat <- read.table("microglial/filter_1e-05_THA_eur_expression_exposure_data.txt",header = T)
head(test_eqtl_exp_dat)

test_eqtl_exp_dat$id.exposure <- test_eqtl_exp_dat$SNP

outcome<-read.table("microglial/overlap_THA_nallsEtAl2019.txt",header = T)
head(outcome)
outcome$samplesize <- outcome$N_cases + outcome$N_controls
write.table(outcome,file="microglial/final_overlap_THA_eQTL_nallsEtAl2019.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

outcome_dat <- read_outcome_data(
  filename = "microglial/final_overlap_THA_eQTL_nallsEtAl2019.txt",sep = "\t",snp_col = "SNP",
  beta_col = "Beta",se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2",
  pval_col = "P",samplesize_col = "samplesize"
)



dat <- harmonise_data(test_eqtl_exp_dat, outcome_dat)


##Heterogeneity test
het <- mr_heterogeneity(dat)
het

##There was no heterogeneity in the results


##Pleiotropy test
pleio <- mr_pleiotropy_test(dat)
pleio



###keep doing MR analysis

res <- mr(dat)
res <- subset_on_method(res)

0.05/296
#0.0001689189
write.table(res,file="Microglial_THA_PD_MR_res.txt", row.names=F, quote=F,
            col.names=T, sep="\t")

##plot

bleeding <- read.table("microglial/anno_Microglial_THA_PD_MR_res.txt", head=TRUE, sep = "\t")
head(bleeding)
library(qqman)
require("RColorBrewer")
display.brewer.all(type = "qual")

manhattan(bleeding,chr="CHROM",bp="POS",p="pval",snp="exposure",
          col =c("#ec4035", "#1799d1"),
          main = "Microglial_PD_MR_res",genomewideline =-log10(0.0001689189),
          suggestiveline = F,annotatePval = 0.0001689189)
