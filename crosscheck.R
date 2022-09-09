library(tidyverse)
library(data.table)
library(readxl)
library(fasttime)

# -1. Helpers ----

drg <- function(name, dc, db) (!is.na(drug_codes[Name==name,Item_Code]) & grepl(drug_codes[Name==name,Item_Code], dc, ignore.case=T)) | (!is.na(drug_codes[Name==name,BNF]) & grepl(drug_codes[Name==name,BNF], db, ignore.case=T))

covar_dz <- function(ct, name, icd) {
  hx1820 <- dx_clean_[grepl(icd, code, ignore.case=T), unique(patient_pssn)]
  hx21 <- dx_latest_index[date < index.date & grepl(icd, code, ignore.case=T), unique(patient_pssn)]
  hx <- unique(c(hx1820, hx21))
  ct[, c(name) := as.numeric(patient_pssn %in% hx)]
  cat(length(hx),"\n")
  return(ct)
}

bl_dz <- function(ct, dz_def) {
  for(i in 1:nrow(dz_def)) {
    cat(dz_def[i,Name], "...")
    ct <- covar_dz(ct, dz_def[i,Name], dz_def[i,Regex])
  }
  return(ct)
}

covar_drug <- function(ct, name, within=90) {
  hx21 <- drug_latest_index[presc_start < index.date & presc_end >= (index.date-within) & drg(name, item_cd, bnfno_p), unique(patient_pssn)]
  hx <- unique(c(hx21))
  coln <- if(within==90) name else paste0(name, ".", within)
  ct[, c(coln) := as.numeric(patient_pssn %in% hx)]
  cat(length(hx),"\n")
  return(ct)
}

cc.match <- function(cohort, case.var, date.var, date.gap, match.var, K, seed, norep_strict=F) {
  setorderv(cohort, c(date.var, match.var, "patient_pssn"))
  cases <- cohort[get(case.var)==1]
  ctrls <- cohort[get(case.var)==0]
  
  matched <- copy(cases[, c("patient_pssn", match.var, date.var), with=F])
  .Random.seed <- seed
  
  c_d <- NA
  pool_case <- NULL
  
  for(i in 1:nrow(cases)) {
    if(i%%1000 == 0) cat(i,"\n")
    case <- cases[i]
    if(is.na(c_d) | case[[date.var]] != c_d) {
      c_d <- case[[date.var]]
      ctrls_ <- ctrls[abs(get(date.var)-c_d) <= date.gap]
    }
    
    if(is.null(pool_case) | !all(unlist(lapply(c(date.var, match.var), function(v) case[[v]]==pool_case[[v]])))) {
      pool_case <- copy(case)
      pool <- Reduce(intersect, lapply(match.var, function(v) {c_v <- case[[v]]; ctrls_[get(v)==c_v, unique(patient_pssn)]}))
    }
    
    if(length(pool)==0) {cat("No match for", case$patient_pssn, nrow(ctrls_), case$Age, case$sex, "\n"); ctrl <- rep(NA, K); print(case$index.date)}
    else ctrl <- c(sample(as.character(pool), min(length(pool),K), replace=F), rep(NA, K-min(length(pool),K)))
    matched[i, paste0("M",1:K) := as.list(as.numeric(ctrl))]
    matched[i, pool.size := length(pool)]
  }
  
  return(matched)
}

output_matches <- function(md, ct, K=10) {
  ms <- melt(cbind(match=1:nrow(md), md[, c("patient_pssn", paste0("M",1:K)), with=F]), id.vars=c("match"), value.name="patient_pssn")
  ct.m <- merge(ct, ms[, .(patient_pssn, match)], by="patient_pssn")
  return(ct.m)
}



# 0. Load data ----

# cohort
cohort <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/4.cohort_full.RDS")
cohort[death_date_ymd=="", death_date_ymd := NA]
cohort[, c(grep("^Age[.]|^Sex[.]", names(cohort), value=T)) := NULL]


# Lab COVID
lab_all_covid <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/LAB_ALL_COVID.RDS")
lab_chp_covid <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/wanning effect/CHP confirmed cases/LAB_CHP_COVID.rds")
lab_chp_covid[, test.date := date][, date := pmin(test.date, report.date, na.rm=T)][, exact_date := as.numeric(!is.na(test.date) & date==test.date)]
lab_chp_covid <- merge(lab_chp_covid, cohort[, .(PseudoID, patient_pssn)], by="PseudoID", all.x=T)


# DX
load("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/DX.RData")
dx_latest[, date := as.Date(fastPOSIXct(date))]

# Attendance latest
load("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/A10X.RData")
admin_latest[, date := as.Date(fastPOSIXct(date))]

# RX
drug_latest <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/RX_latest.RDS")
drug_latest[, presc_start := as.Date(fastPOSIXct(presc_start_date_ymd))]
drug_latest[, presc_end := presc_start + presc_duration_day]



# 1. Identify cases ----

# All cases within study period ---
# CHP confirmed case
cases_chp <- lab_chp_covid %>% filter(test_pos=="PCR" & date >= as.Date("2022-01-01") & date <= as.Date("2022-03-31")) %>% 
  arrange(date, test_pos)# 727021
cases_chp <- cases_chp %>% group_by(PseudoID) %>% slice(1)  # 726552
cases_chp <- cases_chp %>% filter(!is.na(patient_pssn)) # 492075
cases_chp <- as.data.table(cases_chp)
# HA PCR +ve
cases_ha <- lab_all_covid[grepl("^21[3567]$",T_NUM) & result=="detected" & date >= as.Date("2022-01-01") & date <= as.Date("2022-03-31")]
setorder(cases_ha, patient_pssn, date, order, -T_NUM)
cases_ha <- cases_ha[, .SD[1], by=patient_pssn]  # 177512
cases_ha <- cases_ha[patient_pssn %in% cohort$patient_pssn]  # 175228

cases_pcr <- rbind(cases_chp[, .(patient_pssn, date, test_type = test_pos, variant, CTl30, exact_date)],
                   cases_ha[, .(patient_pssn, date, test_type = paste0("HA_",T_NUM))],
                   fill=T) 
setorder(cases_pcr, patient_pssn, date, test_type, exact_date)
cases_pcr <- cases_pcr[, .SD[1], by=patient_pssn]  # 514000

# RAT
cases_rat <- lab_chp_covid %>% filter(test_pos=="RAT" & date >= as.Date("2022-01-01") & date <= as.Date("2022-03-31")) %>% 
  arrange(PseudoID, -exact_date, date, test_pos) %>% 
  group_by(PseudoID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(!is.na(patient_pssn)) %>% 
  filter(!patient_pssn %in% cases_pcr$patient_pssn) # 231832

cases_full <- rbind(cases_pcr, select(cases_rat, patient_pssn, date, test_type=test_pos, variant, CTl30, exact_date))  # 745832

# Death and vaccination status
cases_full <- merge(cases_full, cohort, by="patient_pssn")
# ---

# Hospitalization after infection
cases_adm <- merge(admin_latest, cases_full[, .(patient_pssn, covid.date=date)], by="patient_pssn")
cases_adm <- cases_adm[Source=="1.IP" & date>=covid.date]
setorder(cases_adm, patient_pssn, date, order)
cases_adm <- cases_adm[, .SD[1], by=patient_pssn]  # 35709
cases_full <- merge(cases_full, cases_adm[, .(patient_pssn, ip.adm_date=as.Date(date), ip.dischg_date=as.Date(dischg_date_ymd))], by="patient_pssn", all.x=T)
rm(cases_adm)
# ---

# ICU admission after infection
cases_icu <- merge(admin_latest, cases_full[, .(patient_pssn, covid.date=date)], by="patient_pssn")
cases_icu <- cases_icu[i_icu=="1" & date>=covid.date]
setorder(cases_icu, patient_pssn, date, order)
cases_icu <- cases_icu[, .SD[1], by=patient_pssn]  # 704
cases_full <- merge(cases_full, cases_icu[, .(patient_pssn, icu.adm_date=as.Date(date), icu.dischg_date=as.Date(dischg_date_ymd))], by="patient_pssn", all.x=T)
rm(cases_icu)
# ---

# Ventilatory support after infection
cases_vent <- merge(dx_latest[grepl("^39.65|^89.18|^93.9[056]|^96.7|^96.04", code)], cases_full[, .(patient_pssn, covid.date=date)], by="patient_pssn") 
cases_vent <- cases_vent[date>=covid.date]  # 1390
setorder(cases_vent, patient_pssn, date, order)
cases_vent <- cases_vent[, .SD[1], by=patient_pssn]  # 1052
cases_full <- merge(cases_full, cases_vent[, .(patient_pssn, vent.date=as.Date(date))], by="patient_pssn", all.x=T)
rm(cases_vent)
# ---

# Complications after infection
cases_comp <- merge(dx_latest, cases_full[, .(patient_pssn, covid.date = date)], by="patient_pssn")
comp_codes <- setDT(read_excel("codes.xlsx", sheet="comp"))[!is.na(Name)]
for(i in 1:nrow(comp_codes)) {
  cat(comp_codes[i,Name], "...")
  cases_comp_i <- cases_comp[grepl(comp_codes[i,Regex], code, ignore.case=T) & date>=covid.date]
  setorder(cases_comp_i, patient_pssn, date, order, Source, ranking)
  cases_comp_i <- cases_comp_i[, .SD[1], by=patient_pssn]
  cat(nrow(cases_comp_i),"\n")
  cases_full <- merge(cases_full, setnames(cases_comp_i[, .(patient_pssn, date)], "date", paste0(comp_codes[i,Name], ".date")), by="patient_pssn", all.x=T)
}
rm(i, cases_comp_i, cases_comp, comp_codes)
# ---

# Outcome definitions
cases_full[, covid.date := date]  # COVID infection
cases_full[as.Date(death_date_ymd) <= (date+28) & as.Date(death_date_ymd) >= date, covid_death.date := as.Date(death_date_ymd)]  # Death
cases_full[ip.adm_date <= (date+28), covid_ip.date := ip.adm_date]  # Hospitalisation
cases_full[icu.adm_date <= (date+28) | vent.date <= (date+28), covid_icu_vent.date := pmin(icu.adm_date, vent.date, na.rm=T)]  # ICU + ventilatory support 
cases_full[comp.mi.date <= (date+28) | comp.stroke_isch.date <= (date+28), covid_mace.date := pmin(comp.mi.date, comp.stroke_isch.date, na.rm=T)]  # MACE
cases_full[comp.cad.date <= (date+28) | comp.stroke.date <= (date+28) | comp.hf.date <= (date+28), covid_cvd.date := pmin(comp.cad.date, comp.stroke.date, comp.hf.date, na.rm=T)]  # CVD
# ---

OUTCOME <- "covid_icu_vent" # covid, covid_ip, covid_death, covid_icu_vent
cases <- cases_full %>% filter(test_type!="RAT") %>% filter(!is.na(get(paste0(OUTCOME, ".date")))) %>% mutate(index.date = get(paste0(OUTCOME, ".date")))

ctrls_full <- admin_latest %>% filter(date >= as.Date("2022-01-01")) %>% 
  arrange(patient_pssn, date, order, Source) %>% 
  group_by(patient_pssn) %>% 
  slice(1) %>% ungroup() %>% 
  filter(patient_pssn %in% cohort$patient_pssn)
ctrls_full <- ctrls_full[, .SD[1], by=patient_pssn]   # 1798490

ctrls <- ctrls_full %>% filter(!patient_pssn %in% filter(cases_full, !is.na(get(paste0(OUTCOME, ".date"))))$patient_pssn)

cohort_cc <- rbind(select(mutate(cases, case = 1), patient_pssn, index.date, case, Source=test_type), 
                   select(mutate(ctrls, case = 0), patient_pssn, index.date=date, case, Source))
cohort_cc <- merge(cohort_cc, cohort, by="patient_pssn")
cohort_cc$case %>% table()
# covid 514000-1418554; 514000-210563; 745832-1418554
# ip     34142-1763626; 34142-210563; 34864-1763626
# icu     1448-1797038; 1448-210563; 1452-1797038
# death   7740-1790662; 7740-210563; 7879-1790662
# cvd     4735-1793292; 4735-210563; 5346-1793292