##生存
library('survival')
tmp_cli <- read.csv(path_to/example/survival/Test_survival.csv,header = T,row.names = 1)
tmp_cli_new <- tmp_cli
tmp_cli_new <- tmp_cli_new[!is.na(tmp_cli_new$age_at_diagnosis),]
tmp_cli_new <- tmp_cli_new[which(tmp_cli_new$OS.time > 90),]
tmp_cli_new$gender_new <- 0
tmp_cli_new[which(tmp_cli_new$gender == "male"),8] <- 1
tmp_cli_new$age_new <- 1
tmp_cli_new[which(tmp_cli_new$age_at_diagnosis > median(tmp_cli_new$age_at_diagnosis)),9] <- 2
tmp_cli_new <- tmp_cli_new[!is.na(tmp_cli_new$tumor_stage),]
unique(tmp_cli_new$tumor_stage)
tmp_cli_new <- tmp_cli_new[-which(tmp_cli_new$tumor_stage == "not reported"),]
tmp_cli_new$tumor_stage_new <- tmp_cli_new$tumor_stage
tmp_cli_new$tumor_stage_new <- factor(tmp_cli_new$tumor_stage_new,levels = c("stage i","stage ii","stage iii","stage iv"))
tmp_cli_new[grep("stage ia",tmp_cli_new$tumor_stage),10] <- "stage i"
tmp_cli_new[grep("stage ib",tmp_cli_new$tumor_stage),10] <- "stage i"
tmp_cli_new[grep("stage iia",tmp_cli_new$tumor_stage),10] <- "stage ii"
tmp_cli_new[grep("stage iib",tmp_cli_new$tumor_stage),10] <- "stage ii"
tmp_cli_new[grep("stage iic",tmp_cli_new$tumor_stage),10] <- "stage ii"
tmp_cli_new[grep("stage iiia",tmp_cli_new$tumor_stage),10] <- "stage iii"
tmp_cli_new[grep("stage iiib",tmp_cli_new$tumor_stage),10] <- "stage iii"
tmp_cli_new[grep("stage iiic",tmp_cli_new$tumor_stage),10] <- "stage iii"
tmp_cli_new[grep("stage iva",tmp_cli_new$tumor_stage),10] <- "stage iv"
tmp_cli_new[grep("stage ivb",tmp_cli_new$tumor_stage),10] <- "stage iv"
tmp_cli_new$tumor_stage_group <- 1
tmp_cli_new[which(tmp_cli_new$tumor_stage_new %in% c("stage iii","stage iv")),11] <-  2
tmp_cli_new<-tmp_cli_new[!duplicated(tmp_cli_new$X_PATIENT),]
rownames(tmp_cli_new) <- tmp_cli_new$X_PATIENT
tmp_exp <- read.csv(path_to/example/survival/Test_PSI.csv,header = T,row.names = 1)
colnames(tmp_exp)<-gsub("_",'-', colnames(tmp_exp))
tmp_exp$gene_name <- paste0(tmp_exp$symbol,"_",tmp_exp$`as-id`)
rownames(tmp_exp) <- tmp_exp$gene_name
tmp_exp <- tmp_exp[,intersect(colnames(tmp_exp),tmp_cli_new$X_PATIENT)]
tmp_cli_new <- tmp_cli_new[intersect(colnames(tmp_exp),tmp_cli_new$X_PATIENT),]
result_multivariable <- c()
 for(j in 1:nrow(tmp_exp)){
    
    
    cox_list_multi <- list(exp = as.numeric(tmp_exp[j,]),
                           follow_time = tmp_cli_new$OS.time,
                           status = tmp_cli_new$OS,
                           #stage = tmp_cli_new$tumor_stage_group,
                           age = tmp_cli_new$age_new,
                           gender = tmp_cli_new$gender_new
    )
    
    cox_multi <- coxph(Surv(follow_time,status) ~ exp+age+gender, cox_list_multi)
    cox_terms_multi <- summary(cox_multi)
    
    cox_terms_multi$coefficients
    
    
    tmp_result_multivariable <- c(rownames(tmp_exp)[j],cox_terms_multi$coefficients[1,2],cox_terms_multi$coefficients[1,5])
    
    result_multivariable <- rbind(result_multivariable,tmp_result_multivariable)
    
    
  }
  
  colnames(result_multivariable) <- c("gene_name","HR","P")
  
  write.table(result_multivariable, file = paste0("path_to/multivariable_cox/",cancer_name[i],"_cox.txt"),sep = "\t",
              quote = F)