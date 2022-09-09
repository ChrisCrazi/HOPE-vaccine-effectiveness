
library(data.table)
library(readxl)
library(writexl)


# 0. Helpers ----

table1 <- function(cohort, strata) {
  cohort_ <- copy(cohort)
  cohort_[, (grep("^(dx|px|ip|rx|hx)",names(cohort_),value=T)) := lapply(.SD, as.logical), .SDcols=grep("^(dx|px|ip|rx|hx)",names(cohort_),value=T)]
  
  tb1_def <- setDT(read_excel("codes.xlsx", sheet="table1"))
  t1 <- CreateTableOne(tb1_def[!is.na(Name), Name], c(strata), data=cohort_)
  print(t1, smd=T, test=F)
  return(t1)
}

as.data.table.TableOne <- function(t) {
  tb1_def <- setDT(read_excel("codes.xlsx", sheet="table1"))[!is.na(Name), .(Name, Description)]
  tb1_def <- rbind(list(Name="n", Description="n"), tb1_def)
  t <- as.data.frame(print(t, test=F, dropEqual=T, noSpaces=T))
  varlabels <- rownames(t)
  t$Name = sub("^([a-zA-Z0-9._]+).*$", "\\1", varlabels)
  t <- merge(as.data.table(t), tb1_def, by="Name", all.x=T, sort=F)
  t$Description = paste0(t$Description, sub("^([a-zA-Z0-9._]+)", "", varlabels))
  t[!is.na(Description), Name:=Description][, Description:=NULL]
  return(t)
}

as.data.table.clogit <- function(r, n=3) {
  t <- cbind(OR=exp(r$coefficients), exp(confint(r)))
  t <- as.data.table(t, keep.rownames=T)
  #t$rn <- gsub("vaccinated.brand", "Vaccinated: ", t$rn)
  t[, `OR (95% CI)` := paste0(format(round(OR,n),nsmall=n,trim=T), " (", format(round(`2.5 %`,n),nsmall=n,trim=T), " - ", format(round(pmin(`97.5 %`,9999),n),nsmall=n,trim=T), ")")]
  return(t[, .(` `=rn, `OR (95% CI)`)])
}

calcVE <- function(r, a, b, n=1) {
  rbindlist(lapply(lapply(r, function(x) x$model[[a]]), function(x) setnames(data.table::transpose(as.data.table((1-c(OR=exp(x$coefficients[[paste0("vacc.",b)]]), exp(confint(x))[paste0("vacc.",b),]))*100, keep.rownames=T), make.names="V1"), c("VE","Upper","Lower"))), idcol="Days")[, .(Days, `VE (95% CI)`=paste0(format(round(VE,n),nsmall=n,trim=T), " (", format(round(pmax(`Lower`,-9999),n),nsmall=n,trim=T), "-", format(round(Upper,n),nsmall=n,trim=T), ")"), VE, Upper, Lower)]
}


# 1. Load data ----

OUTCOME <- "covid_icu_vent" # covid, covid_ip, covid_death, covid_icu_vent
# SA_NAME <- ""  # -SA_

cohort.matched.bl <- readRDS(paste0("5.cohort_matched.", OUTCOME, "-1.RDS"))
cohort.matched.bl <- check
# cohort.matched.bl <- readRDS(paste0("5.cohort_matched.", OUTCOME, "-SA_RAT.RDS"))
cohort.matched.bl[, `:=`(`Date of vaccination.1st`=as.Date(`Date of vaccination.1st`), `Date of vaccination.2nd`=as.Date(`Date of vaccination.2nd`), `Date of vaccination.3rd`=as.Date(`Date of vaccination.3rd`), `Date of vaccination.4th`=as.Date(`Date of vaccination.4th`))]

# count only on/before event
cohort.matched.bl[`Date of vaccination.1st` > index.date, `:=`(`Date of vaccination.1st`=NA, `Vaccine Brand.1st`=NA)]
cohort.matched.bl[`Date of vaccination.2nd` > index.date, `:=`(`Date of vaccination.2nd`=NA, `Vaccine Brand.2nd`=NA)]
cohort.matched.bl[`Date of vaccination.3rd` > index.date, `:=`(`Date of vaccination.3rd`=NA, `Vaccine Brand.3rd`=NA)]

#cohort.matched.bl[`Vaccine Brand.1st`!=`Vaccine Brand.2nd`, .N]  # =0
cohort.matched.bl[, vacc.1x_biontech := as.numeric(!is.na(`Date of vaccination.1st`) & is.na(`Date of vaccination.2nd`) & `Vaccine Brand.1st`=="BioNTech/Fosun")]
cohort.matched.bl[, vacc.1x_sinovac := as.numeric(!is.na(`Date of vaccination.1st`) & is.na(`Date of vaccination.2nd`) & `Vaccine Brand.1st`=="Sinovac")]
cohort.matched.bl[, vacc.2x_biontech := as.numeric(!is.na(`Date of vaccination.2nd`) & is.na(`Date of vaccination.3rd`)  & `Vaccine Brand.2nd`=="BioNTech/Fosun")]
cohort.matched.bl[, vacc.2x_sinovac := as.numeric(!is.na(`Date of vaccination.2nd`) & is.na(`Date of vaccination.3rd`) & `Vaccine Brand.2nd`=="Sinovac")]
cohort.matched.bl[, vacc.3x_biontech := as.numeric(!is.na(`Date of vaccination.3rd`) & `Vaccine Brand.3rd`=="BioNTech/Fosun" & `Vaccine Brand.2nd`=="BioNTech/Fosun")]
cohort.matched.bl[, vacc.3x_sinovac := as.numeric(!is.na(`Date of vaccination.3rd`) & `Vaccine Brand.3rd`=="Sinovac" & `Vaccine Brand.2nd`=="Sinovac")]
cohort.matched.bl[, vacc.3x_bbs := as.numeric(!is.na(`Date of vaccination.3rd`) & `Vaccine Brand.3rd`=="Sinovac" & `Vaccine Brand.2nd`=="BioNTech/Fosun")]
cohort.matched.bl[, vacc.3x_ssb := as.numeric(!is.na(`Date of vaccination.3rd`) & `Vaccine Brand.3rd`=="BioNTech/Fosun" & `Vaccine Brand.2nd`=="Sinovac")]
# cohort.matched.bl[, vacc.3rd_biontech := as.numeric(!is.na(`Date of vaccination.3rd`) & `Vaccine Brand.3rd`=="BioNTech/Fosun")]
# cohort.matched.bl[, vacc.3rd_sinovac := as.numeric(!is.na(`Date of vaccination.3rd`) & `Vaccine Brand.3rd`=="Sinovac")]


# SA: if <14 days from vaccine dose, count only previous dose
# cohort.matched.bl[(vacc.3x_biontech==1 | vacc.3x_sinovac==1 | vacc.3x_bbs==1 | vacc.3x_ssb==1) & `Date of vaccination.3rd` > (index.date-14), `:=`(vacc.2x_biontech=as.numeric(vacc.3x_biontech|vacc.3x_bbs), vacc.2x_sinovac=as.numeric(vacc.3x_sinovac|vacc.3x_ssb), vacc.3x_biontech=0, vacc.3x_sinovac=0, vacc.3x_bbs=0, vacc.3x_ssb=0)]
# cohort.matched.bl[(vacc.2x_biontech==1 | vacc.2x_sinovac==1) & `Date of vaccination.2nd` > (index.date-14), `:=`(vacc.1x_biontech=vacc.2x_biontech, vacc.1x_sinovac=vacc.2x_sinovac, vacc.2x_biontech=0, vacc.2x_sinovac=0)]
# cohort.matched.bl[(vacc.1x_biontech==1 | vacc.1x_sinovac==1) & `Date of vaccination.1st` > (index.date-14), `:=`(vacc.1x_biontech=0, vacc.1x_sinovac=0)]
# 
# SA: limit vaccine exposure: not count if vaccine date > 180 days from index
# cohort.matched.bl[(vacc.1x_biontech==1 | vacc.1x_sinovac==1) & `Date of vaccination.1st` <= (index.date-180), `:=`(vacc.1x_biontech=0, vacc.1x_sinovac=0)]
# cohort.matched.bl[(vacc.2x_biontech==1 | vacc.2x_sinovac==1) & `Date of vaccination.2nd` <= (index.date-180), `:=`(vacc.2x_biontech=0, vacc.2x_sinovac=0)]
# cohort.matched.bl[(vacc.3x_biontech==1 | vacc.3x_sinovac==1 | vacc.3x_bbs==1 | vacc.3x_ssb==1) & `Date of vaccination.3rd` <= (index.date-180), `:=`(vacc.3x_biontech=0, vacc.3x_sinovac=0, vacc.3x_bbs=0, vacc.3x_ssb=0)]



# factor vaccination status
table(cohort.matched.bl[, rowSums(.SD), .SDcols=grep("^vacc[.]", names(cohort.matched.bl), value=T)], useNA="always")
cohort.matched.bl$id <- seq_along(cohort.matched.bl$match)
cohort.matched.bl <- merge(cohort.matched.bl, unique(melt(cohort.matched.bl[, c("id", "match", grep("^vacc[.]", names(cohort.matched.bl), value=T)), with=F], c("id","match"))[value==1, .(id, match, vacc_status=sub("vacc.","",variable))]), by=c("id","match"), all.x=T)
table(cohort.matched.bl$vacc_status, useNA="always")

# Duration since nth dose of vaccine
cohort.matched.bl[, time.since.1st_dose := as.numeric(index.date-`Date of vaccination.1st`)]
cohort.matched.bl[, time.since.2nd_dose := as.numeric(index.date-`Date of vaccination.2nd`)]
cohort.matched.bl[, time.since.3rd_dose := as.numeric(index.date-`Date of vaccination.3rd`)]
cohort.matched.bl[, time.since.recent_dose := as.numeric(index.date - pmax(`Date of vaccination.1st`, `Date of vaccination.2nd`, `Date of vaccination.3rd`, na.rm=T))]
# if(grepl("-14", SA_NAME)) cohort.matched.bl[, time.since.recent_dose := ifelse(vacc.1x_biontech==1|vacc.1x_sinovac==1, time.since.1st_dose, ifelse(vacc.2x_biontech==1|vacc.2x_sinovac==1, time.since.2nd_dose, time.since.3rd_dose))]  # use only with - SA: if <14 days from vaccine dose, count only previous dose
cohort.matched.bl[, .(min=min(time.since.recent_dose), max=max(time.since.recent_dose)), keyby=vacc_status]

# Days since study start
cohort.matched.bl[, calendar_days := as.numeric(index.date - min(index.date))]



# 2. Descriptives ----
library(tableone)

get_tb1 <- function(ct) {
  ts <- list(table1(ct, c("case")), table1(ct, c("vacc_status")), table1(ct[case==1], c("vacc_status")), table1(ct[case==0], c("vacc_status")))
  names(ts) <- paste0(c("case-ctrl","vacc_status","case-vacc","ctrl-vacc"), " (", deparse(substitute(ct)), ")")
  return(sapply(ts, as.data.table, simplify=F))
}

get_tb2 <- function(ct) {
  ts <- list(ct[, .(Case=sum(case==1), Control=sum(case==0)), keyby=vacc_status])
  names(ts) <- paste0(c("events"), " (", deparse(substitute(ct)), ")")
  return(ts)
}



# 3. Conditional logit ----
library(survival)
#data <- copy(cohort.matched.bl)

#f0 <- "case ~ vacc.1x_biontech + vacc.1x_sinovac + vacc.2x_biontech + vacc.2x_sinovac + vacc.3x_biontech + vacc.3x_sinovac + vacc.3x_bbs + vacc.3x_ssb"
f0_ <- c("1x_biontech", "1x_sinovac", "2x_biontech", "2x_sinovac", "3x_biontech", "3x_sinovac", "3x_bbs", "3x_ssb")

adj <- c("dx.cancer", "dx.crf", "dx.respdz", "dx.dm", "dx.cvd", "dx.dementia", "rx.ras", "rx.bb", "rx.ccb", "rx.diuretic", "rx.nitrates", "rx.lipid", "rx.insulin", "rx.dm", "rx.oac", "rx.apt", "rx.immunsupp")
adj1 <- c("dx.cancer", "dx.crf", "dx.respdz", "dx.dm", "dx.cvd", "dx.dementia", "rx.ras", "rx.bb", "rx.ccb", "rx.diuretic", "rx.nitrates", "rx.lipid", "rx.insulin", "rx.dm", "rx.oac", "rx.apt", "rx.immunsupp", "rx.antibiotics.7", "rx.antiviral.7")

clean_adj <- function(adj, data) {
  data.table(var=adj, N=unlist(lapply(adj, function(a) data[get(a)==1, .N])))[N>0, var]
}


get_res <- function(f,a,d) {
  res <- list(clogit(as.formula(paste(f, "+ strata(match)")), d),
              clogit(as.formula(paste(f, "+", paste(a, collapse=" + "), "+ strata(match)")), d))
  names(res) <- c(paste0(deparse(substitute(f)), " (", deparse(substitute(d)), ")"),
                  paste0(deparse(substitute(f)), ".", deparse(substitute(a)), " (", deparse(substitute(d)), ")"))
  res.tab <- sapply(res, as.data.table, simplify=F)
  return(list(model=res, table=res.tab))
}

#r0 <- get_res(f0, adj, data)
#r1 <- get_res(f0, adj1, data)

library(writexl)
#write_xlsx(c(get_tb1(data), get_tb2(data), r0$table, r1$table), paste0("res.", OUTCOME, SA_NAME, ".xlsx"))

#save(r0, r1, file=paste0("res.", OUTCOME, SA_NAME, ".RData"))



# Timing since vaccination

# if 2-dose-only
data <- copy(cohort.matched.bl[(vacc.2x_biontech==1|vacc.2x_sinovac==1) | is.na(vacc_status)])

# data <- data[Age>=18 & Age<65]
data <- data[Age>=65]

evt.timing <- list()
res.timing <- list()
for(p in list(c(0,13), c(14,30), c(31,60), c(61,90), c(91,120), c(121,150), c(151,180), c(180,Inf))) {
  cat(paste0(p[1],"-",p[2]), "...\n")
  dt <- copy(data[is.na(vacc_status) | (time.since.recent_dose >= p[1] & time.since.recent_dose <= p[2])])
  dt <- dt[match %in% dt[case==1, match]]
  dt <- dt[!match %in% dt[, sum(case==0), by=match][V1==0, match]]
  ft <- paste("case ~", paste(paste0("vacc.", intersect(f0_, dt[,.N, keyby=vacc_status][N>0, vacc_status])), collapse=" + "))
  adj_ <- clean_adj(adj, dt)
  adj1_ <- clean_adj(adj1, dt)
  r0 <- get_res(ft, adj_, dt)
  r1 <- get_res(ft, adj1_, dt)
  res.timing[[paste0(p[1],"-",p[2])]] <- list(model=c(r0$model, r1$model), table=c(r0$table, r1$table))
  evt.timing[[paste0(p[1],"-",p[2])]] <- get_tb2(dt)
}
rm(p, dt, ft, r0, r1)

merge_ <- function(L, by=NULL) Reduce(function(x,y) merge(x, y, all=T, sort=F), L)
res.timing.summary <- c(
                        list(events = merge_(lapply(names(res.timing), function(n) setnames(evt.timing[[n]]$`events (dt)`[c(NA,f0_), .(vacc_status, paste(Case,"/",Control))], c("vacc_status", paste(n,"days")))), by="vacc_status")),
                        lapply(list("crude"="ft (dt)", "adj"="ft.adj_ (dt)", "adj1"="ft.adj1_ (dt)"), function(x) merge_(lapply(names(res.timing), function(n) setnames(res.timing[[n]]$table[[x]], c("vacc_status", paste(n,"days")))), by="vacc_status")) )

tab.timing <- sapply(list("tab.adj"="adj", "tab.adj1"="adj1"), function(a) {
    t <- lapply(list("2x_biontech", "2x_sinovac"), 
              function(n) merge_(list(data.table::transpose(res.timing.summary$events, keep.names="Days", make.names="vacc_status")[, .(Days=sub(" days","",Days), Case=tstrsplit(get(n)," / ")[[1]], Control=tstrsplit(get(n)," / ")[[2]])], 
                                      data.table::transpose(res.timing.summary$crude, keep.names="Days", make.names="vacc_status")[, .(Days=sub(" days","",Days), `Crude OR (95% CI)`=gsub(" - ","-",get(paste0("vacc.",n))))], 
                                      data.table::transpose(res.timing.summary[[a]], keep.names="Days", make.names="vacc_status")[, .(Days=sub(" days","",Days), `Adjusted OR (95% CI)`=gsub(" - ","-",get(paste0("vacc.",n))))],
                                      calcVE(res.timing, paste0("ft.",a,"_ (dt)"), n))))
    rbind(c("BNT162b2", rep(list(""),ncol(t[[1]])-1)), t[[1]], c("CoronaVac", rep(list(""),ncol(t[[2]])-1)), t[[2]], fill=T)
  }, USE.NAMES = T, simplify = F)

# write_xlsx(c(res.timing.summary, tab.timing), paste0("res.timing.", OUTCOME, SA_NAME, ".xlsx"))



# Age subgroups

#r0.l18 <- get_res(f0_, clean_adj(adj, data[Age<18]), data[Age<18])  # only 1 case dead
# r0.1264 <- get_res(f0, adj, data[Age>=12 & Age<65])
# r0.ge65 <- get_res(f0, adj, data[Age>=65])

# write_xlsx_ <- function(sheets, filename) {
#   names(sheets) <- gsub("\\[", "(", gsub("\\]", ")", names(sheets)))
#   write_xlsx(sheets, filename)
# }
# #write_xlsx_(c(get_tb1(data[Age<18]), get_tb2(data[Age<18]), r0.l18$table), paste0("res.covid.l18", SA_NAME, ".xlsx"))
# write_xlsx_(c(get_tb1(data[Age>=12 & Age<65]), get_tb2(data[Age>=12 & Age<65]), r0.1264$table), paste0("res.", OUTCOME, ".1264", SA_NAME, ".xlsx"))
# write_xlsx_(c(get_tb1(data[Age>=65]), get_tb2(data[Age>=65]), r0.ge65$table), paste0("res.", OUTCOME, ".ge65", SA_NAME, ".xlsx"))



# # Age adjusted VE
# 
# evt.age <- list()
# res.age <- list()
# for(p in list(c(12,17), c(18,50), c(51,64), c(65,Inf), c(65,79), c(80,Inf))) {
#   cat(paste0(p[1],"-",p[2]), "...\n")
#   dt <- copy(data[Age >= p[1] & Age <= p[2]])
#   dt <- dt[match %in% dt[case==1, match]]
#   dt <- dt[!match %in% dt[, sum(case==0), by=match][V1==0, match]]
#   ft <- paste("case ~", paste(paste0("vacc.", intersect(f0_, dt[,.N, keyby=vacc_status][N>0, vacc_status])), collapse=" + "))
#   adj_ <- clean_adj(adj, dt)
#   adj1_ <- clean_adj(adj1, dt)
#   r0 <- get_res(ft, adj_, dt)
#   r1 <- get_res(ft, adj1_, dt)
#   res.age[[paste0(p[1],"-",p[2])]] <- list(model=c(r0$model, r1$model), table=c(r0$table, r1$table))
#   evt.age[[paste0(p[1],"-",p[2])]] <- get_tb2(dt)
# }
# rm(p, dt, ft, r0, r1)
# 
# merge_ <- function(L, by=NULL) Reduce(function(x,y) merge(x, y, all=T, sort=F), L)
# res.age.summary <- c(list(events = merge_(lapply(names(res.age), function(n) setnames(evt.age[[n]]$`events (dt)`[c(NA,f0_), .(vacc_status, paste(Case,"/",Control))], c("vacc_status", paste(n,"yo")))), by="vacc_status")),
#                      lapply(lapply(list("crude"="ft (dt)", "adj"="ft.adj_ (dt)", "adj1"="ft.adj1_ (dt)"), function(x) merge_(lapply(names(res.age), function(n) setnames(res.age[[n]]$table[[x]], c("vacc_status", paste(n,"yo")))), by="vacc_status")), function(x) rbind(x[vacc_status %in% paste0("vacc.",f0_)], x[!vacc_status %in% paste0("vacc.",f0_)])) )

# 
# library(metafor)  # Pooled OR using meta-analysis methods with reference population as weights
# pop_full <- setDT(read_excel("ref_pop_2021.xlsx", sheet="full_pop"))[type==1]
# 
# calc_asOR <- function(ra, pop, adj, n=3) {
#   est <- rbindlist(lapply(names(ra), function(n) as.data.table(coef(summary(ra[[n]]$model[[adj]]))[intersect(names(ra[[n]]$model[[adj]]$coefficients), paste0("vacc.",f0_)), c(1,3)], keep.rownames=T)[, age.grp := n]))
#   est <- merge(est, pop[, .(age.grp, age.prop)], by="age.grp")
#   est.rma <- sapply(f0_, function(x) rma(yi=coef, sei=`se(coef)`, weights=age.prop, data=est[rn==paste0("vacc.",x)], method="FE"), simplify=F, USE.NAMES=T)
#   rbindlist(lapply(est.rma, function(e) list(asOR_full=paste0(format(round(exp(e$b),n),nsmall=n), " (", format(round(exp(e$ci.lb),n),nsmall=n), " - ", format(round(pmin(exp(e$ci.ub),9999),n),nsmall=n,trim=T), ")"))), idcol="vacc_status")
# } 
# 
# res.age.asOR <- lapply(list(crude="ft (dt)",adj="ft.adj (dt)",adj1="ft.adj1 (dt)"), function(a) calc_asOR(res.age, pop_full, a))


# write_xlsx(c(res.age.summary), paste0("res.age.", OUTCOME, SA_NAME, ".xlsx"))
# save(evt.age, res.age, res.age.summary, res.age.asOR, file=paste0("res.age.", OUTCOME, SA_NAME, ".RData"))


