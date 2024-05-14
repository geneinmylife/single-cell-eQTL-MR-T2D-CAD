#package
require(data.table)
require(TwoSampleMR)
require(ieugwasr)
require(fs)

#select SNPs from the original file
# List files
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_steiger",
                             pattern = "[.]csv",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_steiger",
                                  pattern = "[.]csv",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

EUR_SNP<-fread("/PATH_TO_MY_FOLDER/1000G_EUR_v3/EUR.bim")

for (i in 1:46){
  dat_SF<-fread(cd4_data_files[i])
  exp_dat<-dat_SF[,c("SNP","chr.exposure","pos.exposure","other_allele.exposure","effect_allele.exposure","pval.exposure","beta.exposure","se.exposure",
                     "eaf.exposure","samplesize.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure")]
  #clumping within each phenotype_id
  exp_dat<-exp_dat[which(exp_dat$SNP%in%EUR_SNP$V2),]
  exp_dat2<-ld_clump(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure), 
                     clump_kb = 10000, 
                     clump_r2 = 0.1,
                     clump_p = 1,
                     pop = "EUR",
                     plink_bin = "/PATH_TO_MY_FOLDER/Plink-1.90/plink",
                     bfile = "/PATH_TO_MY_FOLDER/1000G_EUR_v3/EUR")
  write.csv(exp_dat2,file=paste("/PATH_TO_MY_FOLDER/exposure_use_clumping/",cd4_data_files_name[i],"_clumping.csv",sep=""),row.names = FALSE)
  rm(dat_SF,exp_dat,exp_dat2)
  gc()
}
