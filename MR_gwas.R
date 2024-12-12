#########1.
rm (list = ls ())
gc()
library(TwoSampleMR)
load("MRGWASDATA.RData")
remove(aries_mqtl,gtex_eqtl,metab_qtls,proteomic_qtls,ao)
exp=gwas_catalog
###读取数据
out=read_outcome_data(
  filename = "finngen_R10_C3_PANCREAS_ADENO_DUCTAL_EXALLC.gz",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval",
  se_col = "sebeta",
)
out$outcome="finngen_R10_C3_PANCREAS_ADENO_DUCTAL_EXALLC"
a=subset(gwas_catalog,pval<=5e-15)
b=a[order(a$pval), ]
table(duplicated(b$Phenotype_simple))
c = b[!duplicated(b$Phenotype_simple),]
name = c$Phenotype_simple[1:200]
remove(a,b,c)
name = b$Phenotype_simple
name=name[!duplicated(name)]
name=name[1:200]
table(duplicated(name))
#x=subset(a,SNP %in% out$SNP) 
gc()
final_df <- data.frame(id.exposure = character(),
                       id.outcome = character(),
                       outcome = character(),
                       exposure = character(),
                       method = character(),
                       nsnp= numeric(),
                       b= numeric(),
                       se= numeric(),
                       pval= numeric())
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

for (i in 1:length(name)) {
  tryCatch({
    data <- subset(gwas_catalog, Phenotype == name[i])
    #data <- data[order(data$pval), ]
    data <- subset(data,data$pval<=5e-08)
    exposure_dat <- format_data(data)
    exposure_dat <- clump_data(exposure_dat)#是否clump
    outcome_dat <- subset(out,SNP %in% exposure_dat$SNP )
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    mr_results <- mr(dat)
    final_df <- rbind(final_df, mr_results)
  }, error = function(e) {
    # 输出错误信息，可以根据需要进行处理
    cat("Error occurred in iteration", i, ":", conditionMessage(e), "\n")
  })
} 

library(openxlsx)
write.xlsx(file = "finngen_R10_C3_PANCREAS_ADENO_DUCTAL_EXALLC.xlsx",x = results)
library(dplyr)
filtered_df <- results %>%
  filter(nsnp >= 5 & method=="Inverse variance weighted" &pval<=0.05)
write.xlsx(file = "finngen_R10_C3_PANCREAS_ADENO_DUCTAL_EXALLC_filtere.xlsx",x = filtered_df)















#########1.
rm (list = ls ())
gc()
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
library(TwoSampleMR)
load("MRGWASDATA.RData")
remove(aries_mqtl,gtex_eqtl,metab_qtls,proteomic_qtls,ao)
exp=gwas_catalog
###读取数据
out=read_outcome_data(
  filename = "finngen_R10_J10_VIRALPNEUMO.gz",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval",
  se_col = "sebeta",
)
out$outcome="finngen_R10_J10_VIRALPNEUMO"
a<-subset(out,abs(beta.outcome)>0.6)
b<-subset(gwas_catalog,SNP %in% a$SNP)
c<-as.data.frame(table(b$Phenotype_simple))
d <- c[order(c$Freq,decreasing = T), ]
name=d$Var1
name=as.character(name)
name=name[1:200]#取p值最小的前200个

remove(a,b,c,d)
gc()
final_df <- data.frame(id.exposure = character(),
                       id.outcome = character(),
                       outcome = character(),
                       exposure = character(),
                       method = character(),
                       nsnp= numeric(),
                       b= numeric(),
                       se= numeric(),
                       pval= numeric())

for (i in 1:length(name)) {
  tryCatch({
    data <- subset(gwas_catalog, Phenotype == name[i])
    #data <- data[order(data$pval), ]
    data <- subset(data,data$pval<=5e-08)
    exposure_dat <- format_data(data)
    exposure_dat <- clump_data(exposure_dat)#是否clump
    outcome_dat <- subset(out,SNP %in% exposure_dat$SNP )
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    mr_results <- mr(dat)
    final_df <- rbind(final_df, mr_results)
  }, error = function(e) {
    # 输出错误信息，可以根据需要进行处理
    cat("Error occurred in iteration", i, ":", conditionMessage(e), "\n")
  })
} 

library(openxlsx)
write.xlsx(file = "finngen_R10_J10_VIRALPNEUMO.xlsx",x = final_df)
library(dplyr)
filtered_df <- final_df %>%
  filter(nsnp >= 5 & method=="Inverse variance weighted" &pval<=0.05)
write.xlsx(file = "finngen_R10_J10_VIRALPNEUMO_filtere.xlsx",x = filtered_df)

##################################
