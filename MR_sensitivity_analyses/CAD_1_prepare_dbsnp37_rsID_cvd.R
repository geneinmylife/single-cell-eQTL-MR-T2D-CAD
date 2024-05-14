#package
require(data.table)
require(fs)
require(stringr)
require(purrr)
require(tidyr)
require(dplyr)

#select SNPs from the original file
# List files
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_rsID",
                             pattern = "[.]csv",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use_rsID",
                                  pattern = "[.]csv",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

#dbSNP data
rsID<-fread("/PATH_TO_MY_FOLDER/37/bg37_rsID.txt")
rsID<-rsID[,c("ID","POS")]
gc()

for (i in 1:length(cd4_data_files)){
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
  #merge with bg37
  cd4_full_use<-merge(cd4_full_use,rsID,by="ID")
  write.csv(cd4_full_use,file=paste("/PATH_TO_MY_FOLDER/exposure_use_cvd_rsID/",cd4_data_files_name[i],"_bg37.csv",sep=""),row.names=FALSE)
  rm(cd4_full,cd4_full_dup,cd4_full_dupix,cd4_full_nodup,cd4_full_use,eaix1,eaix2)
  gc()
}


