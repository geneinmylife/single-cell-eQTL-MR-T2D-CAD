#package
require(data.table)
require(TwoSampleMR)
#require(ieugwasr)
require(fs)

#select SNPs from the original file
# List files
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_rsID",
                             pattern = "[.]csv",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_rsID",
                                  pattern = "[.]csv",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

#data
t2d<-fread("/PATH_TO_MY_FOLDER/outcome/DIAMANTE-EUR.sumstat.txt")

for (i in 1:46){
  cd4_full<-fread(cd4_data_files[i])
  cd4_full<-cd4_full[which(!is.na(cd4_full$ID)),]
  #check duplicate
  cd4_full$ChrPosPheno<-paste(cd4_full$chr,cd4_full$position,cd4_full$phenotype_id,sep=":")
  cd4_full_dupix<-cd4_full[which(duplicated(cd4_full$ChrPosPheno)),]
  cd4_full_nodup<-cd4_full[-which(cd4_full$ChrPosPheno%in%cd4_full_dupix$ChrPosPheno),]
  cd4_full_dup<-cd4_full[which(cd4_full$ChrPosPheno%in%cd4_full_dupix$ChrPosPheno),]
  cd4_full_dup$check<-0
  eaix1<-which(cd4_full_dup$A1==cd4_full_dup$REF & cd4_full_dup$A2==cd4_full_dup$ALT)
  eaix2<-which(cd4_full_dup$A1==cd4_full_dup$ALT & cd4_full_dup$A2==cd4_full_dup$REF)
  cd4_full_dup$check[eaix1]<-1
  cd4_full_dup$check[eaix2]<-1
  cd4_full_dup<-cd4_full_dup[which(cd4_full_dup$check==1),]
  cd4_full_dup<-cd4_full_dup[,-"check"]
  cd4_full_use<-rbind(cd4_full_nodup,cd4_full_dup)
  #merge with outcome GWAS
  cd4_full_use<-merge(cd4_full_use,t2d,by.x="ID",by.y="rsID")
  cd4_full_use$effect_allele<-toupper(cd4_full_use$effect_allele)
  cd4_full_use$other_allele<-toupper(cd4_full_use$other_allele)
  #exposure EAF
  cd4_full_use$eaf_exp<-cd4_full_use$maf
  eafix1<-which(cd4_full_use$A2==cd4_full_use$effect_allele & cd4_full_use$effect_allele_frequency>0.5)
  eafix2<-which(cd4_full_use$A1==cd4_full_use$effect_allele & cd4_full_use$effect_allele_frequency<0.5)
  cd4_full_use$eaf_exp[eafix1]<-1-cd4_full_use$maf[eafix1]
  cd4_full_use$eaf_exp[eafix2]<-1-cd4_full_use$maf[eafix2]
  cd4_full_use$sample_exp<-119
  cd4_full_use$ncase_out<-74124
  cd4_full_use$ncontrol_out<-824006
  exp_dat<-format_data(cd4_full_use,
                       type = "exposure",
                       phenotype_col = "phenotype_id",
                       snp_col = "ID",
                       beta_col = "slope",
                       se_col = "slope_se",
                       eaf_col = "eaf_exp",
                       effect_allele_col = "A2",
                       other_allele_col = "A1",
                       pval_col = "pval_nominal",
                       samplesize_col = "sample_exp",
                       chr_col = "chr",
                       pos_col = "position")
  out_dat<-format_data(cd4_full_use,
                       type = "outcome",
                       phenotype_col = "phenotype_id",
                       snp_col = "ID",
                       beta_col = "Fixed-effects_beta",
                       se_col = "Fixed-effects_SE",
                       eaf_col = "effect_allele_frequency",
                       effect_allele_col = "effect_allele",
                       other_allele_col = "other_allele",
                       pval_col = "Fixed-effects_p-value",
                       ncase_col = "ncase_out",
                       ncontrol_col = "ncontrol_out",
                       chr_col = "chr",
                       pos_col = "position")
  dat<-harmonise_data(exposure_dat = exp_dat,outcome_dat = out_dat,action = 2)
  dat<-dat[which(dat$mr_keep==TRUE & dat$outcome==dat$exposure),]
  dat_SF<-steiger_filtering(dat)
  dat_SF<-dat_SF[which(dat_SF$steiger_dir==TRUE),]
  write.csv(dat_SF,file=paste("/PATH_TO_MY_FOLDER/exposure_use_steiger/",cd4_data_files_name[i],"_steiger.csv",sep=""),row.names=FALSE)
  #exp_dat<-dat_SF[,c("SNP","chr.exposure","pos.exposure","other_allele.exposure","effect_allele.exposure","pval.exposure","beta.exposure","se.exposure",
  #                   "eaf.exposure","samplesize.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure")]
  #clumping within each phenotype_id
  #exp_dat<-clump_data(exp_dat,clump_r2 = 0.1)
  #write.csv(exp_dat,file=paste("/PATH_TO_MY_FOLDER/exposure_use_clumping/",cd4_data_files_name[i],"_clumping.csv",sep=""),row.names=FALSE)
  rm(cd4_full,cd4_full_dup,cd4_full_dupix,cd4_full_nodup,cd4_full_use,exp_dat,out_dat,dat,dat_SF,eaix1,eaix2,eafix1,eafix2)
  gc()
}
