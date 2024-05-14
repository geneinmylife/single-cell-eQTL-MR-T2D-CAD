#package
require(data.table)
require(TwoSampleMR)
#require(ieugwasr)
require(fs)

#select SNPs from the original file
# List files
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_cvd_rsID",
                             pattern = "[.]csv",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_cvd_rsID",
                                  pattern = "[.]csv",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

#data
cvd<-fread("/PATH_TO_MY_FOLDER/outcome/CAD_GWAS_primary_discovery_meta.tsv")
cvd<-cvd[,c("CHR","BP","Allele1","Allele2","Freq1","Effect","StdErr","P-value","Cases","N")]
cvd$Allele1<-toupper(cvd$Allele1)
cvd$Allele2<-toupper(cvd$Allele2)
colnames(cvd)<-c("chr","POS","EA","NEA","EAF","Effect_cvd","SE_cvd","P_cvd","Cases_cvd","N_cvd")
gc()

for (i in 1:length(cd4_data_files)){
  cd4_full<-fread(cd4_data_files[i])
  cd4_full<-cd4_full[which(!is.na(cd4_full$ID)),]
  #merge with outcome GWAS
  cd4_full<-merge(cd4_full,cvd,by=c("chr","POS"))
  #check duplicate
  #cd4_full$ChrPosPheno<-paste(cd4_full$chr,cd4_full$position,cd4_full$phenotype_id,sep=":")
  cd4_full_dupix<-cd4_full[which(duplicated(cd4_full$ChrPosPheno)),]
  cd4_full_nodup<-cd4_full[-which(cd4_full$ChrPosPheno%in%cd4_full_dupix$ChrPosPheno),]
  cd4_full_dup<-cd4_full[which(cd4_full$ChrPosPheno%in%cd4_full_dupix$ChrPosPheno),]
  cd4_full_dup$check<-0
  eaix1<-which(cd4_full_dup$A1==cd4_full_dup$EA & cd4_full_dup$A2==cd4_full_dup$NEA)
  eaix2<-which(cd4_full_dup$A1==cd4_full_dup$NEA & cd4_full_dup$A2==cd4_full_dup$EA)
  cd4_full_dup$check[eaix1]<-1
  cd4_full_dup$check[eaix2]<-1
  cd4_full_dup<-cd4_full_dup[which(cd4_full_dup$check==1),]
  cd4_full_dup<-cd4_full_dup[,-"check"]
  cd4_full_use<-rbind(cd4_full_nodup,cd4_full_dup)
  #exposure EAF
  cd4_full_use$eaf_exp<-cd4_full_use$maf
  eafix1<-which(cd4_full_use$A2==cd4_full_use$EA & cd4_full_use$EAF>0.5)
  eafix2<-which(cd4_full_use$A1==cd4_full_use$EA & cd4_full_use$EAF<0.5)
  cd4_full_use$eaf_exp[eafix1]<-1-cd4_full_use$maf[eafix1]
  cd4_full_use$eaf_exp[eafix2]<-1-cd4_full_use$maf[eafix2]
  cd4_full_use$sample_exp<-119
  cd4_full_use$ncase_out<-cd4_full_use$Cases_cvd
  cd4_full_use$ncontrol_out<-cd4_full_use$N_cvd-cd4_full_use$Cases_cvd
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
                       pos_col = "POS")
  out_dat<-format_data(cd4_full_use,
                       type = "outcome",
                       phenotype_col = "phenotype_id",
                       snp_col = "ID",
                       beta_col = "Effect_cvd",
                       se_col = "SE_cvd",
                       eaf_col = "EAF",
                       effect_allele_col = "EA",
                       other_allele_col = "NEA",
                       pval_col = "P_cvd",
                       ncase_col = "ncase_out",
                       ncontrol_col = "ncontrol_out",
                       chr_col = "chr",
                       pos_col = "POS")
  dat<-harmonise_data(exposure_dat = exp_dat,outcome_dat = out_dat,action = 2)
  dat<-dat[which(dat$mr_keep==TRUE & dat$outcome==dat$exposure),]
  dat_SF<-steiger_filtering(dat)
  dat_SF<-dat_SF[which(dat_SF$steiger_dir==TRUE),]
  #exp_dat<-dat_SF[,c("SNP","chr.exposure","pos.exposure","other_allele.exposure","effect_allele.exposure","pval.exposure","beta.exposure","se.exposure",
  #                   "eaf.exposure","samplesize.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure")]
  #clumping within each phenotype_id
  #exp_dat<-clump_data(exp_dat,clump_r2 = 0.1)
  write.csv(dat_SF,file=paste("/PATH_TO_MY_FOLDER/exposure_use_cvd_steiger/",cd4_data_files_name[i],"_rsID_steiger.csv",sep=""),row.names=FALSE)
  #write.csv(exp_dat,file=paste("/PATH_TO_MY_FOLDER/exposure_use_clumping/",cd4_data_files_name[i],"_rsID_clumping.csv",sep=""),row.names=FALSE)
  rm(cd4_full,cd4_full_dup,cd4_full_dupix,cd4_full_nodup,cd4_full_use,exp_dat,out_dat,dat,dat_SF,eaix1,eaix2,eafix1,eafix2)
  gc()
}
