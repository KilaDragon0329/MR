#########单独处理某个暴露########
rm (list = ls ())
gc()
library(TwoSampleMR)
load("MRGWASDATA.RData")

remove(aries_mqtl,gtex_eqtl,metab_qtls,proteomic_qtls,ao)
exp_id <- read.csv("lbr/exp_id.txt", header=F, sep="\t")
exp_id=exp_id[,1]

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
#out$outcome="finngen_R10_AB1_MALARIA_NOS"
for (i in 1:length(exp_id)) {
  exp=subset(gwas_catalog,gwas_catalog$Phenotype_simple==exp_id[i])
  exp <- exp[order(exp$pval), ]
  #exp =exp[]
  exp <- head(exp, 20)
  #exp=exp[exp$pval>5e-08,]
  exposure_dat <- format_data(exp)
  outcome_dat <- subset(out,SNP %in% exposure_dat$SNP )
  dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)  
  mrResult <- mr(dat)
  #final_df <- rbind(final_df, mrResult)
  #write.csv(dat,file = "lbr/dat.csv",row.names = F)
  #write.csv(mrResult, file="lbr/table.MRresult.csv", row.names=F)
  mrTab=generate_odds_ratios(mrResult)#OR
  #write.csv(mrTab, file="lbr/table.MRresult_or.csv", row.names=F)
  heterTab=mr_heterogeneity(dat)#异质性
  #write.csv(heterTab, file="lbr/table.heterogeneity.csv", row.names=F)
  pleioTab=mr_pleiotropy_test(dat)#多效性检验
  #write.csv(pleioTab, file="lbr/table.pleiotropy.csv", row.names=F)
  
  png(file=paste0("lbr/pic.scatter_plot_",exp_id[i],".png"), width=480, height=520)  #绘制散点图
  mr_scatter_plot(mrResult, dat)
  dev.off()
  
  res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
  png(file=paste0("lbr/pic.forest_",exp_id[i],".png"), width=480, height=520)#pic.forest
  mr_forest_plot(res_single)
  dev.off()
  
  png(file=paste0("lbr/ic.funnel_plot_",exp_id[i],".png"), width=480, height=520)#pic.funnel_plot
  mr_funnel_plot(singlesnp_results = res_single)
  dev.off()
  
  png(file=paste0("lbr/pic.leaveoneout_",exp_id[i],".png"), width=480, height=520)#pic.leaveoneout
  mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
  dev.off()
  }



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
    data <- subset(exp, Phenotype == name[i])
    data <- data[order(data$pval), ]
    data <- head(data, 20)
    exposure_dat <- format_data(data)
    #exposure_dat <- clump_data(exposure_dat)#是否clump
    outcome_dat <- subset(out,SNP %in% exposure_dat$SNP )
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    mr_results <- mr(dat)
    final_df <- rbind(final_df, mr_results)
  }, error = function(e) {
    # 输出错误信息，可以根据需要进行处理
    cat("Error occurred in iteration", i, ":", conditionMessage(e), "\n")
  })
}

write.csv(dat,file = "lbr/dat.csv",row.names = F)
#孟德尔随机化分析
mrResult=mr(dat,)
write.csv(mrResult, file="lbr/table.MRresult.csv", row.names=F)

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="lbr/table.MRresult_or.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="lbr/table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="lbr/table.pleiotropy.csv", row.names=F)

#绘制散点图
png(file="lbr/pic.scatter_plot.png", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
png(file="lbr/pic.forest.png", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
png(file="lbr/pic.funnel_plot.png", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
png(file="lbr/pic.leaveoneout.png", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()
