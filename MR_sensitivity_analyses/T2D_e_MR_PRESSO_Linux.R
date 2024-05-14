#package
require(data.table)
require(arrow)
require(fs)

require(stringr)
require(purrr)
require(tidyr)
require(dplyr)

require(MRPRESSO)

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
for (i in 1:46){
  cd4_data<-fread(cd4_data_files[i])
  cd4_IV<-read.csv(paste("/PATH_TO_MY_FOLDER/exposure_use_clumping/",cd4_data_files_name[i],"_clumping.csv",sep=""))
  cd4_data2<-merge(cd4_data,cd4_IV,by.x=c("SNP","id.exposure"),by.y=c("rsid","id"))
  cd4_data2<-cd4_data2[,-43]
  #2 or more SNPs
  cd4_data_check<-cd4_data2[which(duplicated(cd4_data2$id.exposure)),]
  cd4_data_unique<-cd4_data2[-which(cd4_data2$id.exposure%in%cd4_data_check$id.exposure),]
  cd4_data_2plus<-cd4_data2[-which(cd4_data2$id.exposure%in%cd4_data_unique$id.exposure),]
  #MR-PRESSO for 2 or more SNPs
  cd4_data_2plus$outcome<-"t2d"
  exposure_names<-unique(cd4_data_2plus$exposure)
  for (j in 1:length(exposure_names)){
    data_for_MR<-as.data.frame(cd4_data_2plus[which(cd4_data_2plus$exposure==exposure_names[j]),])
    if(nrow(data_for_MR)>=4){
      # Run MR-PRESSO global method
      MR_PRESSO_results<-mr_presso(BetaOutcome = "beta.outcome", 
                                   BetaExposure = "beta.exposure", 
                                   SdOutcome = "se.outcome", 
                                   SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE, 
                                   DISTORTIONtest = TRUE, 
                                   data = data_for_MR, 
                                   NbDistribution = 1000,  
                                   SignifThreshold = 0.05)
      MR_PRESSO_results<-cbind(exposure_names[j],nrow(data_for_MR),
                               MR_PRESSO_results$`Main MR results`[1,2],MR_PRESSO_results$`Main MR results`[1,3],MR_PRESSO_results$`Main MR results`[1,4],MR_PRESSO_results$`Main MR results`[1,6],
                               MR_PRESSO_results$`Main MR results`[2,2],MR_PRESSO_results$`Main MR results`[2,3],MR_PRESSO_results$`Main MR results`[2,4],MR_PRESSO_results$`Main MR results`[2,6])
      write.table(MR_PRESSO_results,
                  paste("/PATH_TO_MY_FOLDER/results/T2D/MR_PRESSO/",cd4_data_files_name[i],"_MRPRESSO.txt",sep=""),
                  col.names=FALSE,append=TRUE,row.names = FALSE,quote=FALSE)
    }else{
      j<-j+1
    }
  }
}
