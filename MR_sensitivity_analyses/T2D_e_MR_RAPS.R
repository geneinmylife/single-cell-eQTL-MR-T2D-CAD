#package
require(data.table)
require(arrow)
require(fs)

require(stringr)
require(purrr)
require(tidyr)
require(dplyr)

require(TwoSampleMR)
require(mr.raps)

#select SNPs from the original file
# List files
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_steiger",
                             pattern = "[.]csv",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_steiger",
                                  pattern = "[.]csv",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

#input file
for (i in 1:length(cd4_data_files_name)){
  cd4_data<-fread(cd4_data_files[i])
  cd4_IV<-read.csv(paste("/PATH_TO_MY_FOLDER/exposure_use_clumping/",cd4_data_files_name[i],"_clumping.csv",sep=""))
  cd4_data2<-merge(cd4_data,cd4_IV,by.x=c("SNP","id.exposure"),by.y=c("rsid","id"))
  cd4_data2<-cd4_data2[,-43]
  #2 or more SNPs
  cd4_data_check<-cd4_data2[which(duplicated(cd4_data2$id.exposure)),]
  cd4_data_unique<-cd4_data2[-which(cd4_data2$id.exposure%in%cd4_data_check$id.exposure),]
  cd4_data_2plus<-cd4_data2[-which(cd4_data2$id.exposure%in%cd4_data_unique$id.exposure),]
  #MR-RAP for 2 or more SNPs
  cd4_data_2plus$outcome<-"t2d"
  exposure_names<-unique(cd4_data_2plus$exposure)
  for (j in 1:length(exposure_names)){
    cd4_data_2plus_use<-cd4_data_2plus[which(cd4_data_2plus$exposure==exposure_names[j]),]
    try(results_raps<-mr.raps.all(b_exp = cd4_data_2plus_use$beta.exposure,
                                  b_out = cd4_data_2plus_use$beta.outcome,
                                  se_exp = cd4_data_2plus_use$se.exposure,
                                  se_out = cd4_data_2plus_use$se.outcome),silent = TRUE)
    MR_RAPS_results<-cbind(exposure_names[j],nrow(cd4_data_2plus_use),
                           results_raps[1,1],results_raps[1,2],results_raps[1,3],results_raps[1,4])
    write.table(MR_RAPS_results,
                paste("/PATH_TO_MY_FOLDER/results/T2D/2+SNPs_RAPS/try/",cd4_data_files_name[i],"_RAPS.txt",sep=""),
                col.names=FALSE,append=TRUE,row.names = FALSE,quote=FALSE)
  }
}

