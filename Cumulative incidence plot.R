# 7. cumulative incidence plot ----
final_data <- read_rds("mm_final_data.RDS")

final_data$mm_status <- as.numeric(final_data$mm_status)
final_data$mm_status[final_data$mm_status >=4] <- ">=4"
final_data$mm_status <- as.factor(final_data$mm_status)

final_data$treat <- paste0(final_data$vaccinated, final_data$`Vaccine Brand.1st`)

colnames(final_data)<-gsub(", ","_",colnames(final_data))

final_data$mm_cancer <- ifelse(final_data$mm_cancer_lymphoma == 1 |
                                 final_data$mm_cancer_metastatic == 1 |
                                 final_data$`mm_cancer_metastatic & non-metastatic (icpc)` == 1 |
                                 final_data$`mm_cancer_non-metastatic` == 1, 1, 0)  

colnames(final_data)[17:42] <- gsub(" ","_",colnames(final_data))[17:42]
colnames(final_data)[17:42] <- gsub("'","_",colnames(final_data))[17:42]

# obtain the first aesi date to calculate censor status and time
aesi <- read_rds("mm_aesi.RDS")

aesi_first <- merge(x = aesi, y = final_data, by = "patient_pssn", all.x = TRUE)
aesi_first <- aesi_first %>%
  filter(ymd(date) - ymd(`Date of vaccination.1st`) <= 28) %>%
  select(patient_pssn, dx, class, date) %>%
  group_by(patient_pssn) %>%
  arrange(as.Date(date), .by_group = TRUE) %>%
  slice(1)  %>%
  rename(aesi_dx = dx, aesi_class = class, aesi_date = date)

final_data <- merge(x = final_data, y = aesi_first, by = "patient_pssn", all.x = TRUE)

# 0 for censored, 1 for event
final_data$censor <- 0
final_data$censor[!is.na(final_data$aesi_dx) & 
                    (is.na(ymd(final_data$`Date of vaccination.2nd`)) |
                       ymd(final_data$`Date of vaccination.2nd`) > ymd(final_data$aesi_date))] <- 1

final_data$cox_death <- ymd(final_data$death_date_ymd)
final_data$cox_death[is.na(final_data$cox_death)] <- "2021-07-31"
final_data$cox_death <- ymd(final_data$cox_death)

final_data$cox_dose2 <- ymd(final_data$`Date of vaccination.2nd`)
final_data$cox_dose2[is.na(final_data$cox_dose2)] <- "2021-07-31"
final_data$cox_dose2 <- ymd(final_data$cox_dose2)

final_data$cox_event <- ymd(final_data$aesi_date)
final_data$cox_event[is.na(final_data$cox_event)] <- "2021-07-31"
final_data$cox_event <- ymd(final_data$cox_event)

final_data$cox_28 <- ymd(final_data$`Date of vaccination.1st`) + 28

final_data$cox_date <- apply(final_data[,c("cox_death","cox_dose2","cox_event","cox_28")], 1, min)
final_data$time <- ymd(final_data$cox_date) - ymd(final_data$`Date of vaccination.1st`) 

colnames(final_data)<-gsub(", ","_",colnames(final_data))
colnames(final_data)[17:45] <- gsub(" ","_",colnames(final_data))[17:45]
colnames(final_data)[17:45] <- gsub("'","_",colnames(final_data))[17:45]

final_data$vac<-0
final_data$vac[final_data$`Vaccine Brand.1st`=='Sinovac']<-1
final_data$vac[final_data$`Vaccine Brand.1st`=='BioNTech/Fosun']<-2

final_data$surv.cat<-1
final_data$surv.cat[final_data$vac==1&final_data$mm_status==1]<-2
final_data$surv.cat[final_data$vac==2&final_data$mm_status==1]<-3
final_data$surv.cat[final_data$vac==0&final_data$mm_status!=1]<-4
final_data$surv.cat[final_data$vac==1&final_data$mm_status!=1]<-5
final_data$surv.cat[final_data$vac==2&final_data$mm_status!=1]<-6

V3.1 <- data.frame(table(final_data$time[final_data$censor==1 & 
                                           final_data$vac ==0 & 
                                           final_data$mm_status == 1])) %>% 
  rename(V1 = Var1) %>% 
  mutate(V1 = as.numeric(V1)-1) %>% 
  complete(V1 = full_seq(0:28, 1), fill = list(Freq = 0)) %>% 
  rename(V2 = V1) %>% 
  rename(case= Freq) %>% 
  mutate(V1 = 1) %>% 
  relocate(V1, .before = V2) %>% 
  mutate(case= cumsum(case)) %>% 
  mutate(n= nrow(final_data[final_data$vac == 0 & final_data$mm_status == 1]))

V3.2 <- data.frame(table(final_data$time[final_data$censor==1 & 
                                           final_data$vac ==1 & 
                                           final_data$mm_status == 1])) %>% 
  rename(V1 = Var1) %>% 
  mutate(V1 = as.numeric(V1)-1) %>% 
  complete(V1 = full_seq(0:28, 1), fill = list(Freq = 0)) %>% 
  rename(V2 = V1) %>% 
  rename(case= Freq) %>% 
  mutate(V1 = 2) %>% 
  relocate(V1, .before = V2) %>% 
  mutate(case= cumsum(case)) %>% 
  mutate(n= nrow(final_data[final_data$vac == 1 & final_data$mm_status == 1]))

V3.3 <- data.frame(table(final_data$time[final_data$censor==1 & 
                                           final_data$vac ==2 & 
                                           final_data$mm_status == 1])) %>% 
  rename(V1 = Var1) %>% 
  mutate(V1 = as.numeric(V1)-1) %>% 
  complete(V1 = full_seq(0:28, 1), fill = list(Freq = 0)) %>% 
  rename(V2 = V1) %>% 
  rename(case= Freq) %>% 
  mutate(V1 = 3) %>% 
  relocate(V1, .before = V2) %>% 
  mutate(case= cumsum(case)) %>% 
  mutate(n= nrow(final_data[final_data$vac == 2 & final_data$mm_status == 1]))

V3.4 <- data.frame(table(final_data$time[final_data$censor==1 & 
                                           final_data$vac ==0 & 
                                           final_data$mm_status != 1])) %>% 
  rename(V1 = Var1) %>% 
  mutate(V1 = as.numeric(V1)-1) %>% 
  complete(V1 = full_seq(0:28, 1), fill = list(Freq = 0)) %>% 
  rename(V2 = V1) %>% 
  rename(case= Freq) %>% 
  mutate(V1 = 4) %>% 
  relocate(V1, .before = V2) %>% 
  mutate(case= cumsum(case)) %>% 
  mutate(n= nrow(final_data[final_data$vac == 0 & final_data$mm_status != 1]))

V3.5 <- data.frame(table(final_data$time[final_data$censor==1 & 
                                           final_data$vac ==1 & 
                                           final_data$mm_status != 1])) %>% 
  rename(V1 = Var1) %>% 
  mutate(V1 = as.numeric(V1)-1) %>% 
  complete(V1 = full_seq(0:28, 1), fill = list(Freq = 0)) %>% 
  rename(V2 = V1) %>% 
  rename(case= Freq) %>% 
  mutate(V1 = 5) %>% 
  relocate(V1, .before = V2) %>% 
  mutate(case= cumsum(case)) %>% 
  mutate(n= nrow(final_data[final_data$vac == 1 & final_data$mm_status != 1]))

V3.6 <- data.frame(table(final_data$time[final_data$censor==1 & 
                                           final_data$vac ==2 & 
                                           final_data$mm_status != 1])) %>% 
  rename(V1 = Var1) %>% 
  mutate(V1 = as.numeric(V1)-1) %>% 
  complete(V1 = full_seq(0:28, 1), fill = list(Freq = 0)) %>% 
  rename(V2 = V1) %>% 
  rename(case= Freq) %>% 
  mutate(V1 = 6) %>% 
  relocate(V1, .before = V2) %>% 
  mutate(case= cumsum(case)) %>% 
  mutate(n= nrow(final_data[final_data$vac == 2 & final_data$mm_status != 1]))

#https://www.health.pa.gov/topics/HealthStatistics/Statistical-Resources/UnderstandingHealthStats/Documents/Confidence_Intervals_for_a_Crude_Rate.pdf
plots <- bind_rows(V3.1, V3.2, V3.3, V3.4, V3.5, V3.6) %>% 
  mutate(V3 = case/n,
         V4 = (case - 1.96*sqrt(case))/n,
         V5 = (case + 1.96*sqrt(case))/n) %>% 
  select(-c(case, n))

plots1<-plots[plots$V1==1,] 
plots2<-plots[plots$V1==2,] 
plots3<-plots[plots$V1==3,] 
plots4<-plots[plots$V1==4,] 
plots5<-plots[plots$V1==5,] 
plots6<-plots[plots$V1==6,] 

tiff("Figure 2.tiff", height = 8, width = 8.5, units = 'in',compression = "lzw", res = 500)
par(mfrow=c(2,2))

par(mar = c(4, 4, 4, 1.5))
plot(1,type="n",xlim=c(0,28),ylim=c(0,0.005),axes=FALSE,xlab='',ylab='')
axis(2,pos=0,at=c(0, 0.001, 0.002, 0.003, 0.004, 0.005),labels=c('0', '0.001', '0.002', '0.003', '0.004', '0.005'),las=1,cex.axis=1)
axis(1,pos=0,at=c(0:28),labels=c(0:28),cex.axis=1)
title(xlab = 'Follow-up time (in days)', mgp = c(1.5, 1, 0))
title(ylab = 'Cumulative incidence of adverse event of special interest', mgp = c(3, 1, 0))

with(plots1,polygon(c(V2,rev(V2)),c(V4,rev(V5)),
                    col = rgb(.2,.8,.2,alpha=0.5), border = FALSE))
with(plots4,polygon(c(V2,rev(V2)),c(V4,rev(V5)),
                    col = rgb(1,.455,.549,alpha=0.5), border = FALSE))

with(plots1,lines(V2,V3,lty=1,col='darkgreen',lwd=2))
with(plots4,lines(V2,V3,lty=1,col='red',lwd=2))

polygon(c(1,1,3,3),c(.00410,.00400,.00400,.00410),col = rgb(.2,.8,.2,alpha=0.5), border = FALSE)
text(3.5,.00410,'Unvaccinated without multimorbidity',cex=0.8,adj=c(0,1))
polygon(c(1,1,3,3),c(.00450,.00440,.00440,.00450),col = rgb(1,.455,.549,alpha=0.5), border = FALSE)
text(3.5,.00450,'Unvaccinated with multimorbidity',cex=0.8,adj=c(0,1))

par(mar = c(4, 4, 4, 1.5))
plot(1,type="n",xlim=c(0,28),ylim=c(0, 0.005),axes=FALSE,xlab='',ylab='')
axis(2,pos=0,at=c(0, 0.001, 0.002, 0.003, 0.004, 0.005),labels=c('0', '0.001', '0.002', '0.003', '0.004', '0.005'),las=1,cex.axis=1)
axis(1,pos=0,at=c(0:28),labels=c(0:28),cex.axis=1)
title(xlab = 'Follow-up time (in days)', mgp = c(1.5, 1, 0))
title(ylab = 'Cumulative incidence of adverse event of special interest', mgp = c(3, 1, 0))

with(plots2,polygon(c(V2,rev(V2)),c(V4,rev(V5)),
                    col = rgb(1,1,0,alpha=0.5), border = FALSE))
with(plots5,polygon(c(V2,rev(V2)),c(V4,rev(V5)),
                    col = rgb(.2,.2,.8,alpha=0.5), border = FALSE))
with(plots2,lines(V2,V3,lty=1,col='yellow2',lwd=2))
with(plots5,lines(V2,V3,lty=1,col='blue',lwd=2))


polygon(c(1,1,3,3),c(.00410,.00400,.00400,.00410),col = rgb(1,1,0,alpha=0.5), border = FALSE)
text(3.5,.00410,'CoronaVac recipients without multimorbidity',cex=0.8,adj=c(0,1))
polygon(c(1,1,3,3),c(.00450,.00440,.00440,.00450),col = rgb(.2,.2,.8,alpha=0.5), border = FALSE)
text(3.5,.00450,'CoronaVac recipients with multimorbidity',cex=0.8,adj=c(0,1))

par(mar = c(4, 4, 4, 1.5))
plot(1,type="n",xlim=c(0,28),ylim=c(0,0.005),axes=FALSE,xlab='',ylab='')
axis(2,pos=0,at=c(0, 0.001, 0.002, 0.003, 0.004, 0.005),labels=c('0', '0.001', '0.002', '0.003', '0.004', '0.005'),las=1,cex.axis=1)
axis(1,pos=0,at=c(0:28),labels=c(0:28),cex.axis=1)
title(xlab = 'Follow-up time (in days)', mgp = c(1.5, 1, 0))
title(ylab = 'Cumulative incidence of adverse event of special interest', mgp = c(3, 1, 0))

with(plots6,polygon(c(V2,rev(V2)),c(V4,rev(V5)),
                    col = rgb(1,.6,0,alpha=0.5), border = FALSE))
with(plots3,polygon(c(V2,rev(V2)),c(V4,rev(V5)),
                    col = rgb(.5,0,.5,alpha=0.5), border = FALSE))
with(plots6,lines(V2,V3,lty=1,col='orange',lwd=2))
with(plots3,lines(V2,V3,lty=1,col='purple',lwd=2))


polygon(c(1,1,3,3),c(.00410,.00400,.00400,.00410),col = rgb(.5,0,.5,alpha=0.5), border = FALSE)
text(3.5,.00410,'Comirnaty recipients without multimorbidity',cex=0.8,adj=c(0,1))
polygon(c(1,1,3,3),c(.00450,.00440,.00440,.00450),col = rgb(1,.6,0,alpha=0.5), border = FALSE)
text(3.5,.00450,'Comirnaty recipients with multimorbidity',cex=0.8,adj=c(0,1))

dev.off()