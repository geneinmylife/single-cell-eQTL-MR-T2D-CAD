#package
require(data.table)
require(arrow)
require(fs)

require(stringr)
require(purrr)
require(tidyr)
require(dplyr)

#select SNPs from the original file
# List files
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_full",
                             pattern = "[.]parquet",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_full",
                                  pattern = "[.]parquet",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

#set up a P-value threshold (i.e. <0.05), F-statistic threshold (i.e. >=5 based on the meeting with Chris on 28 Feb 2024)
for (i in 1:length(cd4_data_files)){
  cd4_full<-read_parquet(cd4_data_files[i])
  cd4_full<-cd4_full[which(cd4_full$pval_nominal<0.05),]
  cd4_full$f<-(cd4_full$slope/cd4_full$slope_se)^2
  cd4_full<-cd4_full[which(cd4_full$f>=5),]
  cd4_full<-separate(cd4_full,colnames(cd4_full)[2],into=c("chr","position","A1","A2"),sep="_")
  write.csv(cd4_full,file=paste("/PATH_TO_MY_FOLDER/exposure_use/",cd4_data_files_name[i],".csv",sep=""),row.names=FALSE)
  rm(cd4_full)
  gc()
}

