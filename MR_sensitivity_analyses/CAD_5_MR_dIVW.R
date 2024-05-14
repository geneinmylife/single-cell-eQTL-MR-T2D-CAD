#package
require(data.table)
require(arrow)
require(fs)

require(stringr)
require(purrr)
require(tidyr)
require(dplyr)

require(MendelianRandomization)
require(TwoSampleMR)

#select SNPs from the original file
# List files
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_cvd_steiger",
                             pattern = "[.]csv",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_cvd_steiger",
                                  pattern = "[.]csv",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

#input file
for (i in 1:length(cd4_data_files_name)){
  cd4_data<-fread(cd4_data_files[i])
  cd4_IV<-read.csv(paste("/PATH_TO_MY_FOLDER/exposure_use_cvd_clumping/",cd4_data_files_name[i],"_clumping.csv",sep=""))
  cd4_data2<-merge(cd4_data,cd4_IV,by.x=c("SNP","id.exposure"),by.y=c("rsid","id"))
  cd4_data2<-cd4_data2[,-43]
  #2 or more SNPs
  cd4_data_check<-cd4_data2[which(duplicated(cd4_data2$id.exposure)),]
  cd4_data_unique<-cd4_data2[-which(cd4_data2$id.exposure%in%cd4_data_check$id.exposure),]
  cd4_data_2plus<-cd4_data2[-which(cd4_data2$id.exposure%in%cd4_data_unique$id.exposure),]
  #run mr_divw for 2 or more SNPs
  cd4_data_2plus$outcome<-"cvd"
  cd4_data_2plus_mr<-dat_to_MRInput(cd4_data_2plus)
  for (j in 1:length(cd4_data_2plus_mr)){
    data_for_MR<-cd4_data_2plus_mr[[j]]
    try(MR_DIVW<-mr_divw(data_for_MR),silent = TRUE)
    MR_DIVW_results<-cbind(names(cd4_data_2plus_mr[j]),MR_DIVW$SNPs,MR_DIVW$Estimate,MR_DIVW$StdError,MR_DIVW$CILower,MR_DIVW$CIUpper,MR_DIVW$Pvalue)
    write.table(MR_DIVW_results,
                paste("/PATH_TO_MY_FOLDER/results/CVD/2+SNPs_MR-DIVW_full/",cd4_data_files_name[i],"_MRDIVW.txt",sep=""),
                col.names=FALSE,append=TRUE,row.names = FALSE,quote=FALSE)
  }
}


