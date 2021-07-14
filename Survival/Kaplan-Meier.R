library('survival')
library('survminer')
Test<-read.csv("path_to/example/survival/Test_PSI.csv",header = T)
Test_S<-read.csv("path_to/example/survival/Test_survival.csv",header = T)
colnames(Test)<-gsub("_",'-',    colnames(Test))
Test_S_filter<-subset(Test_S,Test_S$OS.time>90)

Test_S_filter_dup<-Test_S_filter[!duplicated(Test_S_filter$X_PATIENT),]
l<- Reduce(intersect, list(colnames(Test),Test_S_filter_dup$X_PATIENT))
Test[,l]->final_Test
cbind(Test[,1:10],final_Test)->final_final_Test
as_event_Test<-as.data.frame(paste(final_final_Test$symbol,final_final_Test$`as-id`,sep="_"))
names(as_event_Test )[1]<-'as_event'
cbind(as_event_Test ,final_final_Test  )  ->final_final_final_Test 
row.names(final_final_final_Test)=final_final_final_Test[,1]
row.names(Test_S_filter_dup)<-Test_S_filter_dup$X_PATIENT
Test_S_final_final<-Test_S_filter_dup[l,]
cbind(Test_S_final_final,t(final_final_final_Test[,-(1:11)]))->all_Test
colnames(all_Test)<-gsub("\\-",'_',colnames(all_Test))


group <- ifelse(all_Test$DPH2_2496>mean(all_Test$DPH2_2496),'high','low')


sfit <- survfit(Surv(OS.time, OS)~group, data=all_Test)
ggsurvplot(sfit, conf.int=F, pval=TRUE,legend.title="PSI level",legend=c(0.7,0.9),legend.labs=c("H_PSI","L_PSI")
           ,risk.table = F,
           title='esample',
           xlab='Time(days)')
