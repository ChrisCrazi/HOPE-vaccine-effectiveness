library(tidyverse)
drug_latest <- readRDS("RDS/RX_latest.RDS")
drug_clean <- readRDS("RDS/RX_clean.RDS")
cohort <- readRDS("RDS/4.cohort_full.RDS")
load("RDS/DX.RData")
LAB_ALL_COVID <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/LAB_ALL_COVID.RDS")
drug_latest %>% filter(str_detect(bnfno_p, "^2\\.12")) %>% .$item_cd %>% table()

pcsk9_users <- drug_latest %>% filter(str_detect(item_cd, "ALIR|EVOL")) %>% .$patient_pssn %>% unique()
statins_users <- drug_latest %>% filter(str_detect(item_cd, "ATOR|ROSU|SIMV")) %>% .$patient_pssn %>% unique()

covid_pts <- LAB_ALL_COVID %>% filter(str_detect(T_NUM, "^21[3567]$") & result=="detected") %>% .$patient_pssn %>% unique()

ascvd_pts <- unique(c(filter(dx_clean, str_detect(code, "^41[0-4]|^43[0-8]|^44[0-3]|^36\\."))$patient_pssn, c(filter(dx_latest, str_detect(code, "^41[0-4]|^43[0-8]|^44[0-3]|^36\\."))$patient_pssn)))
hf_pts <- unique(c(filter(dx_clean, str_detect(code, "^428"))$patient_pssn, c(filter(dx_latest, str_detect(code, "^41[0-4]|^43[0-8]|^44[0-3]|^36\\."))$patient_pssn)))

statins_users %>% length()
pcsk9_users %>% length()
intersect(intersect(ascvd_pts, covid_pts), statins_users) %>% length()
intersect(intersect(ascvd_pts, covid_pts), pcsk9_users) %>% length()
intersect(covid_pts, statins_users) %>% length()
intersect(covid_pts, pcsk9_users) %>% length()
intersect(hf_pts, pcsk9_users) %>% length()

