#package
require(data.table)
require(arrow)
require(fs)

require(stringr)
require(purrr)
require(tidyr)
require(dplyr)

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
for (i in 1:46){
  cd4_data<-fread(cd4_data_files[i])
  cd4_IV<-read.csv(paste("/PATH_TO_MY_FOLDER/exposure_use_cvd_clumping/",cd4_data_files_name[i],"_clumping.csv",sep=""))
  cd4_data2<-merge(cd4_data,cd4_IV,by.x=c("SNP","id.exposure"),by.y=c("rsid","id"))
  cd4_data2<-cd4_data2[,-43]
  #Wald ratio
  cd4_data_check<-cd4_data2[which(duplicated(cd4_data2$id.exposure)),]
  cd4_data_unique<-cd4_data2[-which(cd4_data2$id.exposure%in%cd4_data_check$id.exposure),]
  cd4_data_2plus<-cd4_data2[-which(cd4_data2$id.exposure%in%cd4_data_unique$id.exposure),]
  #run weighted median & weighted mode for 2+ SNPs
  results_wm<-mr(cd4_data_2plus,method_list = "mr_weighted_mode")
  results_wm$outcome<-"cvd"
  write.csv(results_wm,paste("/PATH_TO_MY_FOLDER/results/CVD/2+SNPs_WeightedMode/",cd4_data_files_name[i],"_mode.csv",sep=""),row.names=FALSE)
}

