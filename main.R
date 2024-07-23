library(dplyr)
library(survival)
library(survminer)

clinical <- read.delim("~/Documents/2023_work/2023_lab_course/ym/survival/clinical.cart.2024-07-22/clinical.tsv")
x <- c("case_id","case_submitter_id",'vital_status','days_to_last_follow_up')
clinical <- clinical[,x]
clinical <- clinical %>% distinct()

connect <- read.delim("~/Documents/2023_work/2023_lab_course/ym/survival/gdc_sample_sheet.2024-07-22.tsv")

clinical_con <- merge(clinical, connect, by.x = "case_submitter_id", by.y = "Case.ID")

value <- {}
for(n in c(1:374)){
  file <- file.path("~/Documents/2023_work/2023_lab_course/ym/survival/gdc_download_20240722_124106.302821/",clinical_con$File.ID[n],clinical_con$File.Name[n])
  gex <- read.delim(file, header = TRUE,skip=1)[-c(1:4),]
  value[n] <- gex$fpkm_unstranded[n]
}

median_value <- median(value)
categories <- ifelse(value < median_value, "Low", "High")

sur_data <- data_frame(clinical_con[,c('vital_status','days_to_last_follow_up')],categories)
sur_data <- sur_data[which(sur_data$days_to_last_follow_up != "'--"),]
sur_data <- sur_data[which(sur_data$vital_status != "Not Reported"),]
sur_data$days_to_last_follow_up <- as.numeric(sur_data$days_to_last_follow_up)
vital_binary <- ifelse(sur_data$vital_status=="Dead" , 1, 0)
sur_data$vital_status <- vital_binary


km_fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ categories, data=sur_data)

ggsurvplot(km_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("High", "Low"), legend.title="TSPAN6",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Kaplan-Meier Curve for TSPAN6", 
           risk.table.height=.15)
