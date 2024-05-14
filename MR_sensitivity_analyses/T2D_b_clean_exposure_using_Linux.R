#package
require(data.table)
require(fs)
require(stringr)
require(purrr)
require(tidyr)
require(dplyr)

#select SNPs from the original file
# List files
cd4_data_files <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use",
                             pattern = "[.]csv",
                             full.names = T)
cd4_data_files_name <- list.files(path = "/PATH_TO_MY_FOLDER/exposure_use",
                                  pattern = "[.]csv",
                                  full.names = F)
cd4_data_files_name <- path_ext_remove(cd4_data_files_name)

#data
rsID<-fread("/user/work/qy17115/MRPREG/38/bg38_rsID.txt")

for (i in 1:length(cd4_data_files)){
  cd4_full<-fread(cd4_data_files[i])
  cd4_full$chr<-as.character(cd4_full$chr)
  cd4_full<-merge(cd4_full,rsID,by.x=c("chr","position"),by.y=c("#CHROM","POS"),all.x=TRUE)
  #cd4_full<-cd4_full[which(!is.na(cd4_full$ID)),]
  #check duplicate
  #cd4_full$ChrPosPheno<-paste(cd4_full$chr,cd4_full$position,cd4_full$phenotype_id,sep=":")
  #cd4_full_dupix<-cd4_full[which(duplicated(cd4_full$ChrPosPheno)),]
  #cd4_full_nodup<-cd4_full[-which(cd4_full$ChrPosPheno%in%cd4_full_dupix$ChrPosPheno),]
  #cd4_full_dup<-cd4_full[which(cd4_full$ChrPosPheno%in%cd4_full_dupix$ChrPosPheno),]
  #cd4_full_dup$check<-0
  #for (i in 1:length(cd4_full_dup$ID)){
  #  if(cd4_full_dup$A1[i]==cd4_full_dup$REF[i] & cd4_full_dup$A2[i]==cd4_full_dup$ALT[i]){cd4_full_dup$check[i]<-1}
  #  if(cd4_full_dup$A1[i]==cd4_full_dup$ALT[i] & cd4_full_dup$A2[i]==cd4_full_dup$REF[i]){cd4_full_dup$check[i]<-1}
  #}
  #cd4_full_dup<-cd4_full_dup[which(cd4_full_dup$check==1),]
  #cd4_full_dup<-cd4_full_dup[,-"check"]
  #cd4_full_use<-rbind(cd4_full_nodup,cd4_full_dup)
  write.csv(cd4_full,file=paste("/PATH_TO_MY_FOLDER/exposure_use_rsID/",cd4_data_files_name[i],"_rsID.csv",sep=""),row.names=FALSE)
  rm(cd4_full)
  gc()
}




