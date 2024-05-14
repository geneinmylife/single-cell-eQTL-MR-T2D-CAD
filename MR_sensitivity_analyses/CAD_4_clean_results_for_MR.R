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
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_cvd_steiger",
                             pattern = "[.]csv",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_cvd_steiger",
                                  pattern = "[.]csv",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

#original IVs (use the list for T2D or CVD, if needed)
IV_original<-read.csv("/PATH_TO_MY_FOLDER/IV_list.csv")
IV_original<-IV_original[,c(1:3)]
IV_original$file_name<-as.character(paste(IV_original$cell_time,"_500kb_combined_rsID_bg37_rsID_steiger",sep=""))


for (i in 1:46){
  cd4_data<-fread(cd4_data_files[i])
  cd4_IV<-read.csv(paste("/PATH_TO_MY_FOLDER/exposure_use_cvd_clumping/",cd4_data_files_name[i],"_clumping.csv",sep=""))
  cd4_data2<-merge(cd4_data,cd4_IV,by.x=c("SNP","id.exposure"),by.y=c("rsid","id"))
  cd4_data2<-cd4_data2[,-43]
  IV_original_use<-IV_original[which(IV_original$file_name==cd4_data_files_name[i]),]
  cd4_data3<-cd4_data2[which(cd4_data2$exposure%in%IV_original_use$gene),]
  #run MR analysis using cd4_data3
  write.csv(cd4_data3,file=paste("/PATH_TO_MY_FOLDER/exposure_use_cvd_MR/",cd4_data_files_name[i],"_clumping_forMR.csv",sep=""),row.names=FALSE)
  rm(cd4_data,cd4_IV,cd4_data2,IV_original_use,cd4_data3)
  gc()
}
