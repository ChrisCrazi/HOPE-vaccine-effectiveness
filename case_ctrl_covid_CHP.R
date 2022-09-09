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
    gc()
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


pp.match <- function(dz_pop, nodz_pop, K, seed) {
  dz_cases <- dz_pop[case==1]
  nodz_cases <- nodz_pop[case==1]
  nodz_cases_matched <- pp.match_(dz_cases, nodz_cases, c("Age","sex"), K, seed)
  cat("\n", nrow(dz_cases), ":", nrow(nodz_cases), "|", length(nodz_cases_matched), "-", nrow(nodz_pop[case==0]))
  rbind(nodz_cases[PseudoID %in% nodz_cases_matched], nodz_pop[case==0])
}


pp.match_ <- function(ct_ref, ct_tar, match.var, K, seed) {
  setorderv(ct_ref, c(match.var, "PseudoID"))
  setorderv(ct_tar, c(match.var, "PseudoID"))
  gps <- unique(ct_ref[, match.var, with=F])
  info <- merge(ct_ref[, .N, keyby=match.var], ct_tar[, .N, keyby=match.var], by=match.var, all=T)[, K := N.y/N.x]; print(summary(info$K))
  
  res <- character()
  .Random.seed <- seed
  
  for(i in 1:nrow(gps)) {
    cat("\n", gps[i,Age], gps[i,sex])
    gp_ref <- ct_ref[Reduce(`&`,lapply(match.var, function(v) ct_ref[[v]]==gps[i][[v]]))]
    gp_tar <- ct_tar[Reduce(`&`,lapply(match.var, function(v) ct_tar[[v]]==gps[i][[v]]))]
    gp_tar_rand <- sample(gp_tar$PseudoID, pmin(K*nrow(gp_ref), nrow(gp_tar)), replace=F)
    res <- c(res, gp_tar_rand)
    cat("->", nrow(gp_ref), ":", length(gp_tar_rand), "|", nrow(gp_tar))
  }
  
  return(res)
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


# 1. Identify cases ----

#lab_chp_covid[, exact_date := as.numeric(!is.na(date))][is.na(date), date := report.date][, report.date := NULL]

# All cases within study period ---
# CHP confirmed case (PCR only)
cases_chp <- lab_chp_covid[test_pos=="PCR" & date >= as.Date("2022-01-01") & date <= as.Date("2022-03-31")]  # 727021
setorder(cases_chp, PseudoID, -exact_date, date, test_pos)
cases_chp <- cases_chp[, .SD[1], PseudoID]  # 726551
cases_chp <- cases_chp[!is.na(patient_pssn)]  # 492075

# HA PCR+
cases_ha <- lab_all_covid[grepl("^21[3567]$",T_NUM) & result=="detected" & date >= as.Date("2022-01-01") & date <= as.Date("2022-03-31")]
setorder(cases_ha, patient_pssn, date, order, -T_NUM)
cases_ha <- cases_ha[, .SD[1], by=patient_pssn]  # 177512
cases_ha <- cases_ha[patient_pssn %in% cohort$patient_pssn]  # 175228

# If duplicated dates, use the earliest PCR+ date within the study period
cases_pcr <- rbind(cases_chp[, .(patient_pssn, date, test_type = test_pos, variant, CTl30, exact_date)],
                   cases_ha[, .(patient_pssn, date, test_type = paste0("HA_",T_NUM))],
                   fill=T) 
setorder(cases_pcr, patient_pssn, date, test_type, exact_date)
cases_pcr <- cases_pcr[, .SD[1], by=patient_pssn]  # 514000

# CHP confirmed case (RAT only)
cases_rat <- lab_chp_covid[test_pos=="RAT" & date >= as.Date("2022-01-01") & date <= as.Date("2022-03-31")]
setorder(cases_rat, PseudoID, -exact_date, date, test_pos)
cases_rat <- cases_rat[, .SD[1], PseudoID]  # 426938
cases_rat <- cases_rat[!is.na(patient_pssn)]  # 252241
# Only consider as RAT case if no PCR+ record
cases_rat <- cases_rat[!patient_pssn %in% cases_pcr$patient_pssn]  # 231832

# ---
cases_full <- rbind(cases_pcr, cases_rat[, .(patient_pssn, date, test_type=test_pos, variant, CTl30, exact_date)])  # 745832
# save(cases_chp, cases_ha, cases_rat, cases_full, file="checkpoints/case_ctrl_covid_CHP/1.cases_full.RData")
rm(cases_chp, cases_ha, cases_pcr, cases_rat)


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

# saveRDS(cases_full, "checkpoints/case_ctrl_covid_CHP/1.cases_full.RDS")


# SELECT outcome ----
OUTCOME <- "covid_icu_vent" # covid, covid_ip, covid_death, covid_icu_vent
PCR_ONLY <- T
# SA_NAME <- "-SA_RAT"
cases <- copy(cases_full[PCR_ONLY==F | test_type!="RAT"])[!is.na(get(paste0(OUTCOME, ".date")))][, index.date := get(paste0(OUTCOME, ".date"))]



# 2. Identify controls ----
ctrls_full <- admin_latest[date >= as.Date("2022-01-01")]   # IF IP only -> ctrls_full <- ctrls_full[Source=="1.IP"]  # 332108
setorder(ctrls_full, patient_pssn, date, order, Source)
ctrls_full <- ctrls_full[, .SD[1], by=patient_pssn]  # 1805668, IP=187003
ctrls_full <- ctrls_full[patient_pssn %in% cohort$patient_pssn]  # 1798490

# saveRDS(ctrls_full, file="checkpoints/case_ctrl_covid_CHP/1.ctrls_full.RDS")
# # ---
# 
# # SA: Test-negative case control (first PCR- within study period, exclude if also had PCR+)
# ctrls_neg <- lab_all_covid[grepl("^21[3567]$",T_NUM) & result=="not detected" & date >= as.Date("2022-01-01") & date <= as.Date("2022-03-31")]
# setorder(ctrls_neg, patient_pssn, date, order, -T_NUM)
# ctrls_neg <- ctrls_neg[, .SD[1], by=patient_pssn]  # 276307
# ctrls_neg <- ctrls_neg[patient_pssn %in% cohort$patient_pssn]  # 274504
# ctrls_neg <- ctrls_neg[!patient_pssn %in% cases_full$patient_pssn]   # 210563
# ctrls_neg[, `:=`(Source=paste0("HA_",T_NUM), dischg_date_ymd="", dischg_gp="", i_icu="")]
# 
# saveRDS(ctrls_neg, file="checkpoints/case_ctrl_covid_CHP/1.ctrls_neg.RDS")
# # ---

# Remove cases (both PCR/RAT) for this outcome
ctrls <- ctrls_full[!patient_pssn %in% cases_full[!is.na(get(paste0(OUTCOME, ".date"))), patient_pssn]]
# save(cases, ctrls, file=paste0("checkpoints/case_ctrl_covid_CHP/2.cases_ctrls.", OUTCOME, SA_NAME, ".RData"))



# 3. Cases and controls ----
cohort_cc <- rbind(cases[, .(patient_pssn, index.date, case=1, Source=test_type, info=paste(variant, CTl30, exact_date, sep="|"))], 
                   ctrls[, .(patient_pssn, index.date=date, case=0, Source, info=paste(dischg_date_ymd, dischg_gp, i_icu, sep="|"))] )
cohort_cc <- merge(cohort_cc, cohort, by="patient_pssn")
cohort_cc$case %>% table()
# covid 514000-1418554; 514000-210563; 745832-1418554
# ip     34142-1763626; 34142-210563; 34864-1763626
# icu     1448-1797038; 1448-210563; 1452-1797038
# death   7740-1790662; 7740-210563; 7879-1790662
# cvd     4735-1793292; 4735-210563; 5346-1793292

# Remove people with any PCR+/RAT+ before study period
hx_covid_chp <- lab_chp_covid[date < as.Date("2022-01-01") & !is.na(patient_pssn), unique(patient_pssn)]
hx_covid_ha <- lab_all_covid[grepl("^21[3567]$",T_NUM) & result=="detected" & date < as.Date("2022-01-01"), unique(patient_pssn)]
cohort_cc[, hx.covid := as.numeric(patient_pssn %in% union(hx_covid_chp, hx_covid_ha))]
cohort_cc <- cohort_cc[hx.covid==0]
rm(hx_covid_chp, hx_covid_ha)
cohort_cc$case %>% table()
# 513230-1415530; 513230-210154; 744865-1415530
#  34094-1760231; 34094-210154; 34816-1760231
#   1448-1793595; 1440-210154; 1452-1793595
#   7738-1787221; 7738-210154; 7877-1787221
#   4734-1789850; 4734-210154; 5345-1789850

# Remove died before index.date
cohort_cc <- cohort_cc[is.na(as.Date(death_date_ymd)) | as.Date(death_date_ymd) >= index.date]
# 513175-1415530; 513175-210146; 744757-1415530
#  34094-1760231; 34094-210146; 34816-1760231
#   1448-1793595; 1440-210146; 1452-1793595
#   7738-1787221; 7738-210146; 7877-1787221
#   4734-1789850; 4734-210146; 5345-1789850

# Remove 4th dose
cohort_cc <- cohort_cc[is.na(`Date of vaccination.4th`)]
cohort_cc$case %>% table()
# 513174-1415366; 513174-210110; 744755-1415366
#  34094-1760066; 34094-210110; 34816-1760066
#   1448-1793430; 1448-210110; 1452-1793430
#   7738-1787056; 7738-210110; 7877-1787056
#   4734-1789685; 4734-210110; 5345-1789685

# saveRDS(cohort_cc, file=paste0("checkpoints/case_ctrl_covid_CHP/3.cohort_cc.", OUTCOME, SA_NAME, ".RDS"))



# 4. Baseline characteristics ----

cohort_cc[, Age_5yb := floor(Age/5)]

# * History of comobidities (2018 onwards) ---
dx_clean_ <- dx_clean[patient_pssn %in% cohort_cc$patient_pssn]
dx_latest_index <- merge(dx_latest, cohort_cc[, .(patient_pssn, index.date)], by="patient_pssn")

dxpx_codes <- setDT(read_excel("codes.xlsx", sheet="codes"))[!is.na(Name)]
cohort_cc <- bl_dz(cohort_cc, dxpx_codes)
cohort_cc[, score.cci := (dx.mi+dx.chf+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+(dx.dm_com0&!dx.dm_com1)+dx.dm_com1*2+dx.crf*2+(dx.liver_mild&!dx.liver_modsev)+dx.liver_modsev*3+dx.ulcers+dx.ra+dx.aids*6+dx.cancer*2+dx.cancer_mets*6)]
cohort_cc[score.cci==0, score.cci.trunc := "0"]
cohort_cc[score.cci==1 | score.cci==2, score.cci.trunc := "1-2"]
cohort_cc[score.cci==3 | score.cci==4, score.cci.trunc := "3-4"]
cohort_cc[score.cci>=5, score.cci.trunc := ">=5"]

# # cancer-related
# cancer1stdx <- rbind(dx_clean_[grepl(dxpx_codes[Name=="dx.im.cancer", Regex], code)][, date := as.Date(paste0(date,"-01"))], 
#                      dx_latest_index[date < index.date & grepl(dxpx_codes[Name=="dx.im.cancer", Regex], code)][, index.date := NULL])
# setorder(cancer1stdx, patient_pssn, date, order, Source, ranking, code, code_ext)
# cancer1stdx <- cancer1stdx[, .SD[1], by=patient_pssn]
# cohort_cc <- merge(cohort_cc, cancer1stdx[, .(patient_pssn, dx.im.cancer.site1st = code)], by="patient_pssn", all.x=T)
# 
# assertthat::assert_that(dx_latest_index[, min(date)] < as.Date(dx_latest_index[, min(index.date)]-365))
# for(n in c("dx.im.cancer_radx", "dx.im.cancer_mets")) {
#   cat(n, "...")
#   hx <- dx_latest_index[date >= (index.date-365) & date < index.date & grepl(dxpx_codes[Name==n, Regex], code, ignore.case=T), unique(patient_pssn)]
#   cohort_cc[, c(paste0(n,".365")) := as.numeric(patient_pssn %in% hx)]
#   cohort_cc[, c(n) := NULL]
#   cat(length(hx),"\n")
# }
# rm(n, hx, cancer1stdx)
# 
# # ---


# * Medication use within 90 days (unless specified) ---
# RX
gc()
drug_latest <- readRDS("C:/Users/LabPC14CSMPR/Desktop/Chris/HOPE/booster dose and COVID-19 mortality in multimorbidity/RDS/RX_latest.RDS")
drug_latest[, presc_start := as.Date(fastPOSIXct(presc_start_date_ymd))]
drug_latest[, presc_end := presc_start + presc_duration_day]
gc()
drug_latest_index <- merge(drug_latest, cohort_cc[, .(patient_pssn, index.date)], by="patient_pssn")
gc()
drug_codes <- setDT(read_excel("codes.xlsx", sheet="drugs"))[!is.na(Name)]
for(i in 1:nrow(drug_codes)) {
  cat(drug_codes[i,Name], "...")
  cohort_cc <- covar_drug(cohort_cc, drug_codes[i,Name])
  gc()
}
rm(i)

# antiinfectives within 7 days
cohort_cc[, `:=`(rx.antibiotics.90=rx.antibiotics, rx.antiviral.90=rx.antiviral, rx.antibiotics=NULL, rx.antiviral=NULL)]
cohort_cc <- covar_drug(cohort_cc, "rx.antibiotics", 7)
cohort_cc <- covar_drug(cohort_cc, "rx.antiviral", 7)

# antineoplastic or immunosuppresants in past 12 months
assertthat::assert_that(drug_latest_index[, min(presc_start)] < as.Date(drug_latest_index[, min(index.date)]-365))
for(n in c("rx.cytotoxic","rx.immunsupp", "rx.hormone_cancer", "rx.dmards")) {
  cat(n, "...")
  cohort_cc <- covar_drug(cohort_cc, n, within=365)
  gc()
}
rm(n)

# ---

# saveRDS(cohort_cc, file=paste0("4.cohort_cc_bl.", OUTCOME, ".RDS"))



# 5. Matching ----

SUBPOP <- ""

# [Optional] subpopulations
# [dm]
# cohort_cc <- cohort_cc[dx.t2dm==1 | rx.dm==1] 


# # [cancer]
# cohort_cc[, dx.im.cancer_active := as.numeric(dx.im.cancer_mets.365|dx.im.cancer_radx.365|rx.cytotoxic.365|rx.immunsupp.365|rx.hormone_cancer.365)]
# cohort_cc <- cohort_cc[dx.im.cancer==1 & dx.im.cancer_active==1]  # active cancer: 6851-29325; 532-36610; 2064-35134
# cohort_cc <- cohort_cc[dx.im.cancer==1 & dx.im.cancer_active==0]  # non-active cancer: 7363-34395; 332-42391; 1242-41459
# cohort_cc <- cohort_cc[dx.im.cancer==1]                           # cancer any: 14214-63720; 864-79001; 3306-76593
# # no cancer cases, match to any cancer cases 10:1 by Age and sex
# dz_pop <- cohort_cc[dx.im.cancer==1]
# nodz_pop <- cohort_cc[dx.im.cancer==0 & dx.im.cancer_active==0]  
# set.seed(475)
# seed <- .Random.seed
# cohort_cc <- pp.match(dz_pop, nodz_pop, 4, seed)  # no-cancer: (1:10) 142126-1328361, (1:2) 1650-1679505, (1:4) 12925-1655598
# rm(dz_pop, nodz_pop, seed)


# # [immunocompromised]
# cohort_cc[, dx.im_any := as.numeric(dx.im.cancer|dx.im.sot|dx.im.hsct|dx.im.primary|dx.im.inflam|dx.im.asplen|dx.im.esrd|dx.im.aids|dx.im.cancer_radx.365|rx.cytotoxic.365|rx.immunsupp.365|rx.dmards.365)]
# cohort_cc <- cohort_cc[dx.im_any==1]  # immunocompromised: 23656-106804; 1811-132568; 5463-128394
# # not immuno cases, matched to immuno cases 10:1 by Age and sex
# dz_pop <- cohort_cc[dx.im_any==1]
# nodz_pop <- cohort_cc[dx.im_any==0]
# set.seed(475)
# seed <- .Random.seed
# cohort_cc <- pp.match(dz_pop, nodz_pop, 4, seed)  # no-immuno: (1:10) 231012-1308562, (1:2) 3111-1654488, (1:4) 17793-1631672
# rm(dz_pop, nodz_pop, seed)


# Remove history of outcome for complications
# cohort_cc <- cohort_cc[dx.cvd_comp==0]
# cvd-t2dm 408-365165


# Match on age, sex, cci, index date +- 3
set.seed(475)
seed <- .Random.seed
ratio <- min(floor(cohort_cc[case==0, .N]/cohort_cc[case==1, .N]), 10)
gc()
matched_cc <- cc.match(cohort_cc, "case", "index.date", 3, c("Age_5yb","sex","score.cci.trunc"), ratio, seed)  
matched_cc <- matched_cc[pool.size > 0]
cohort.matched <- output_matches(matched_cc, cohort_cc, K=ratio)
cohort.matched$case %>% table()
# 513073-1026002; 512961-1024925; 744608-1489026
# 34061-338381; 33880-330020; 34779-345495
# 1448-14386; 1443-14016; 1452-14426
# 7719-76236; 7604-71668; 7858-77595
# save(cohort_cc, matched_cc, seed, cohort.matched, file=paste0("checkpoints/case_ctrl_covid_CHP/5.cohort_matched.", OUTCOME, SUBPOP, SA_NAME, ".RData"))
# saveRDS(cohort.matched, file=paste0("checkpoints/case_ctrl_covid_CHP/5.cohort_matched.", OUTCOME, SUBPOP, SA_NAME, ".RDS"))


# 6. Export ----
cohort.matched.export <- copy(cohort.matched)
# ... generate new id here
# Remove unneeded cols
cohort.matched.export[, `:=`(patient_pssn=NULL, PseudoID=NULL, pseudo_hkid=NULL, pseudo_other_doc_no=NULL)]
#cohort.matched.export[, c(grep("^Age[.]|^Sex[.]", names(cohort.matched.export), value=T)) := NULL]

# saveRDS(cohort.matched.export, file=paste0("5.cohort_matched.", OUTCOME, ".RDS"))
saveRDS(cohort.matched.export, file=paste0("5.cohort_matched.", OUTCOME, "-1.RDS"))
# saveRDS(cohort.matched.export, file=paste0("5.cohort_matched.", OUTCOME, SA_NAME, ".RDS"))






# ---
# 4. Baseline characteristics 

# aesi_codes <- setDT(read_excel("codes.xlsx", sheet="aesi"))[!is.na(Name)]
# cohort.matched <- bl_aesi(cohort.matched, aesi_codes)
# # Fix
# cohort.matched[, hx.t1dm := NULL]
# cohort.matched <- covar_dz(cohort.matched, "hx.t1dm", aesi("t1dm"))

# admin_all <- rbind(admin_clean[, .(patient_pssn, date=paste0(date,"-01"), order, Source)], admin_latest[, .(patient_pssn, date, order, Source)])
# admin_all[, date := as.Date(fastPOSIXct(date))]
# admin_all <- merge(admin_all, cohort.matched[, .(patient_pssn, index.date)], by="patient_pssn")
# admin_all[, index.date := as.Date(index.date)]
# 
# admin_hx <- admin_all[date < index.date, .(attn.N=.N), by=.(patient_pssn, Source)]
# admin_hx <- dcast(admin_hx, patient_pssn ~ Source, value.var="attn.N", fill=0)
# names(admin_hx) <- c("patient_pssn", paste0("attn.N.",tolower(substr(names(admin_hx)[-1],3,4))))
# admin_hx[, attn.N := attn.N.ip + attn.N.ae + attn.N.go + attn.N.so]
# cohort.matched <- merge(cohort.matched, admin_hx, by="patient_pssn", all.x=T)
# cohort.matched[is.na(attn.N), `:=`(attn.N=0, attn.N.ip=0, attn.N.ae=0, attn.N.go=0, attn.N.so=0)]
# 
# admin_hx <- admin_all[date < index.date, .(attn.date=max(date)), by=.(patient_pssn, Source)]
# admin_hx <- dcast(admin_hx, patient_pssn ~ Source, value.var="attn.date", fill=NA)
# names(admin_hx) <- c("patient_pssn", paste0("attn.date.",tolower(substr(names(admin_hx)[-1],3,4))))
# cohort.matched <- merge(cohort.matched, admin_hx, by="patient_pssn", all.x=T)
# 
# admin_latest_index <- merge(admin_latest, cohort.matched[, .(patient_pssn, index.date)], by="patient_pssn")
# admin_latest_index[, `:=`(date=as.Date(fastPOSIXct(date)), dischg_date_ymd=as.Date(fastPOSIXct(dischg_date_ymd)))]
# ip_28 <- admin_latest_index[Source=="1.IP" & date < index.date & (is.na(dischg_date_ymd) | dischg_date_ymd >= (index.date-28)), unique(patient_pssn)]
# cohort.matched[, ip.28 := as.numeric(patient_pssn %in% ip_28)]

# stroke/ramsay hunt post index
# stroke_postindex <- dx_latest_index[grepl("^43[012]|^433.[012389]1|^43[46]|^437.[01]|^435", code) & (Source=="1.IP" & ranking=="P") & date >= index.date, unique(patient_pssn)]
# ramsay_postindex <- dx_latest_index[grepl("^053.11", code) & (Source=="1.IP" & ranking=="P") & date >= index.date, unique(patient_pssn)]
# cohort.matched[, stroke.postindex := as.numeric(patient_pssn %in% stroke_postindex)]
# cohort.matched[, ramsay.postindex := as.numeric(patient_pssn %in% ramsay_postindex)]


# unused
# merge_periods <- function(t, gap) {
#   cummax.Date <- function(d) as.Date(cummax(as.integer(d)), "1970-01-01")
#   setorder(t, start, end)
#   t[, period.gap := start-cummax(shift(end, fill=start[1]-1))-1]
#   t[, period.group := cumsum(as.numeric(period.gap > gap))]
#   t[, .(period.start = min(start), period.end = max(end)), by=period.group][,-1]
# }
# ip_hx_latest <- admin_latest[Source=="1.IP", .(patient_pssn, start=as.Date(date), end=as.Date(dischg_date_ymd))][, merge_periods(.SD, 0), by=patient_pssn]



# SA: IP only
# cases_adm <- merge(admin_latest, cases_full[, .(patient_pssn, covid.date=date)], by="patient_pssn")
# cases_adm <- cases_adm[Source=="1.IP" & date>=(covid.date-2) & date<=(covid.date+2)]
# cases_full <- cases_full[patient_pssn %in% cases_adm$patient_pssn] # 24219
# rm(cases_adm)
# cases_full <- cases_full[test_type!="RAT"]  # 23649
# saveRDS(cases_full, "checkpoints/case_ctrl_covid_CHP/1.cases_full-SA_IPneg.RDS")

# SA: Test-negative within 2 days before index (use with IP only)
# ctrls_neg <-  merge(lab_all_covid, ctrls_full[, .(patient_pssn, admin.date = date)], by="patient_pssn")
# ctrls_neg <- ctrls_neg[grepl("^21[3567]$",T_NUM) & result=="not detected" & date>=(admin.date-2) & date<=(admin.date+2)]
# ctrls_full <- ctrls_full[patient_pssn %in% ctrls_neg$patient_pssn]  # 133389
# rm(ctrls_neg)
# saveRDS(ctrls_full, file="checkpoints/case_ctrl_covid_CHP/1.ctrls_full-SA_IPneg.RDS")
# ---

