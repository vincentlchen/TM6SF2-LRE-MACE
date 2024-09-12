library(tidyverse)
library(survminer)
library(survival)
library(tidycmprsk)

setwd("/nfs/corenfs/INTMED-Speliotes-home/vichen/TM6SF2_cirr_CVD")

gen.summ = function(reg, pred, rows) {
  exp = summary(reg)$conf.int[, "exp(coef)"]
  ll = summary(reg)$conf.int[,"2.5%"]
  ul = summary(reg)$conf.int[,"97.5%"]
  p = summary(reg)$coef[,"p-value"]
  summ = paste0(sprintf(exp, fmt='%#.2f'), 
                " (",
                sprintf(ll, fmt='%#.2f'),
                "-",
                sprintf(ul, fmt='%#.2f'),
                ")")
  return(data.frame(pred = pred, summ = summ, p = p)[1:rows,])
}



hr.table = function(t, rowname = "", time1 = "time.mace", status1 = "status.mace", time2 = "time.lre", status2 = "status.lre", cov) {
  
  if (cov == "") {
    cov.tot = "tm6sf2.1 + tm6sf2.2"
  } else {
    cov.tot = paste("tm6sf2.1 + tm6sf2.2", cov, sep = "+")
  }
  cov.all = base::strsplit(cov.tot, "[ ]*\\+[ ]*")[[1]]
  
  # recode status variable to ensure consistent with tidycmprsk
  t[ , status1] = as.factor(case_when(t[, status1] == 1 ~ 2,
                                      t[, status1] == 2 ~ 3,
                                      t[, status1] == 3 ~ 1))
  
  t[ , status2] = as.factor(case_when(t[, status2] == 1 ~ 2,
                                      t[, status2] == 2 ~ 3,
                                      t[, status2] == 3 ~ 1))
  
  # convert from days to years
  t[ , time1] = t[, time1]/365.25
  t[ , time2] = t[, time2]/365.25
  
  # run models
  a = crr(as.formula(paste0("Surv(", time1, ",", status1, ") ~ ", cov.tot)), data = t, failcode = 2)
  b = crr(as.formula(paste0("Surv(", time2, ",", status2, ") ~ ", cov.tot)), data = t, failcode = 2)
  
  # extract data from models
  a1 = as.data.frame(a$tidy)
  b1 = as.data.frame(b$tidy)
  
  table = data.frame(
    name = rep(rowname, 2),
    
    HR1 = c(paste0(sprintf(exp(a1[1,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(a1[1,2] - 1.96*a1[1,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(a1[1,2] + 1.96*a1[1,3]), fmt = "%#.2f"), ")"),
            paste0(sprintf(exp(a1[2,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(a1[2,2] - 1.96*a1[2,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(a1[2,2] + 1.96*a1[2,3]), fmt = "%#.2f"), ")")),
    p1 = signif(c(a1[1,5], a1[2,5]), 3),
    
    HR2 = c(paste0(sprintf(exp(b1[1,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(b1[1,2] - 1.96*b1[1,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(b1[1,2] + 1.96*b1[1,3]), fmt = "%#.2f"), ")"),
            paste0(sprintf(exp(b1[2,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(b1[2,2] - 1.96*b1[2,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(b1[2,2] + 1.96*b1[2,3]), fmt = "%#.2f"), ")")),
    p2 = signif(c(b1[1,5], b1[2,5]), 2),
    n = dim(t[rowSums(is.na(t[, cov.all]))==0, ])[1],
    cases.1 = dim(t[which(rowSums(is.na(t[, cov.all]))==0 & t[, status1] == 2), ])[1],
    cases.2 = dim(t[which(rowSums(is.na(t[, cov.all]))==0 & t[, status2] == 2), ])[1])
  
  hr1=exp(a1[2,2])
  ll1=exp(a1[2,2]-1.96*a1[2,3])
  ul1=exp(a1[2,2]+1.96*a1[2,3])
  
  hr2=exp(b1[2,2])
  ll2=exp(b1[2,2]-1.96*b1[2,3])
  ul2=exp(b1[2,2]+1.96*b1[2,3])
  
  e1 = ifelse((hr1 > 1 & ll1 < 1) | # if beta positive but LL negative -> E=1
                (hr1 < 1 & ul1 > 1) | # if beta negative but UL positive -> E=1
                hr1 == 1, 1,  # if HR == 1 -> E=1
              ifelse(hr1 > 1, 
                     ll1 + sqrt(ll1 * (ll1 - 1)),
                     1/ul1 + sqrt(1/ul1 * (1/ul1 - 1))))
  
  e2 = ifelse((hr2 > 1 & ll2 < 1) | # if beta positive but LL negative -> E=1
                (hr2 < 1 & ul2 > 1) | # if beta negative but UL positive -> E=1
                hr2 == 1, 1,  # if HR == 1 -> E=1
              ifelse(hr2 > 1, 
                     ll2 + sqrt(ll2 * (ll2 - 1)),
                     1/ul2 + sqrt(1/ul2 * (1/ul2 - 1))))
  
  table = mutate(table,
                 E1 = e1,
                 E2 = e2)
  
  return(table)
  
}








cuminc10y = function(data, rowname = "", 
                     pred1 = "tm6sf2", outcome1 = "time.mace", status1 = "status.mace", 
                     pred2 = "tm6sf2", outcome2 = "time.lre", status2 = "status.lre") {
  
  p = cmprsk::cuminc(ftime = data[, outcome1],
                     fstatus = data[, status1],
                     group = data[, pred1])
  
  p0 = as.data.frame(p$"0 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p0 = p0[p0$time1 == max(p0$time1),]
  p0 = p0[dim(p0)[1],]
  p0 = mutate(p0, summ = paste0(sprintf(100*est, fmt = "%#.1f"), "% (",
                                sprintf(100*est^(exp(-1.96*sqrt(var)/(est*log(est)))), fmt = "%#.1f"), "%-",
                                sprintf(100*est^(exp(1.96*sqrt(var)/(est*log(est)))), fmt = "%#.1f"), "%)"))
  
  p1 = as.data.frame(p$"1 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p1 = p1[p1$time1 == max(p1$time1),]
  p1 = p1[dim(p1)[1],]
  p1 = mutate(p1, summ = paste0(sprintf(100*est, fmt = "%#.1f"), "% (",
                                sprintf(100*est^(exp(-1.96*sqrt(var)/(est*log(est)))), fmt = "%#.1f"), "%-",
                                sprintf(100*est^(exp(1.96*sqrt(var)/(est*log(est)))), fmt = "%#.1f"), "%)"))
  
  p2 = as.data.frame(p$"2 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p2 = p2[p2$time1 == max(p2$time1),]
  p2 = p2[dim(p2)[1],]
  p2 = mutate(p2, summ = paste0(sprintf(100*est, fmt = "%#.1f"), "% (",
                                sprintf(100*est^(exp(-1.96*sqrt(var)/(est*log(est)))), fmt = "%#.1f"), "%-",
                                sprintf(100*est^(exp(1.96*sqrt(var)/(est*log(est)))), fmt = "%#.1f"), "%)"))
  
  
  pcomb = NULL
  pcomb$var = rowname
  pcomb = as.data.frame(pcomb)
  pcomb$inc.10y.0.1 = p0$summ
  pcomb$inc.10y.1.1 = p1$summ
  pcomb$inc.10y.2.1 = p2$summ
  pcomb$inc.10y.1.diff = paste0(sprintf(100*(p2$est - p0$est), fmt = "%#.1f"), "%")
  pcomb$events.1 = dim(data[which(data[, outcome1] <= 10*365.25 & data[, status1] == 1), ])[1]
  pcomb$fu.1 = sum(data[data[, outcome1] <= 10*365.25, outcome1]/365.25)
  
  # p value is defined by the CC vs. TT comparison, not the dosage
  data$pred1.mod = ifelse(data[,pred1] == 0, 0, ifelse(data[,pred1] == 2, 1, NA))
  p = cmprsk::cuminc(ftime = data[, outcome1],
                     fstatus = data[, status1],
                     group = data[, "pred1.mod"])
  pcomb$p.cuminc1.mod = ifelse(p$Tests[1,2] < 0.0001, "<0.0001", signif(p$Tests[1,2], 2))
  
  #
  
  p = cmprsk::cuminc(ftime = data[, outcome2],
                     fstatus = data[, status2],
                     group = data[, pred2])
  
  p0 = as.data.frame(p$"0 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p0 = p0[p0$time1 == max(p0$time1),]
  p0 = p0[dim(p0)[1],]
  p0 = mutate(p0, summ = paste0(sprintf(100*est, fmt = "%#.2f"), "% (",
                                sprintf(100*est^(exp(-1.96*sqrt(var)/(est*log(est)))), fmt = "%#.2f"), "%-",
                                sprintf(100*est^(exp(1.96*sqrt(var)/(est*log(est)))), fmt = "%#.2f"), "%)"))
  
  p1 = as.data.frame(p$"1 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p1 = p1[p1$time1 == max(p1$time1),]
  p1 = p1[dim(p1)[1],]
  p1 = mutate(p1, summ = paste0(sprintf(100*est, fmt = "%#.2f"), "% (",
                                sprintf(100*est^(exp(-1.96*sqrt(var)/(est*log(est)))), fmt = "%#.2f"), "%-",
                                sprintf(100*est^(exp(1.96*sqrt(var)/(est*log(est)))), fmt = "%#.2f"), "%)"))
  
  p2 = as.data.frame(p$"2 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p2 = p2[p2$time1 == max(p2$time1),]
  p2 = p2[dim(p2)[1],]
  p2 = mutate(p2, summ = paste0(sprintf(100*est, fmt = "%#.2f"), "% (",
                                sprintf(100*est^(exp(-1.96*sqrt(var)/(est*log(est)))), fmt = "%#.2f"), "%-",
                                sprintf(100*est^(exp(1.96*sqrt(var)/(est*log(est)))), fmt = "%#.2f"), "%)"))
  
  #
  
  pcomb$inc.10y.0.2 = p0$summ
  pcomb$inc.10y.1.2 = p1$summ
  pcomb$inc.10y.2.2 = p2$summ
  pcomb$inc.10y.2.diff = paste0(sprintf(100*(p2$est - p0$est), fmt = "%#.2f"), "%")
  pcomb$p.cuminc2 = ifelse(p$Tests[1,2] < 0.0001, "<0.0001", signif(p$Tests[1,2], 2))
  pcomb$events.2 = dim(data[which(data[, outcome2] <= 10*365.25 & data[, status2] == 1), ])[1]
  pcomb$fu.2 = sum(data[data[, outcome2] <= 10*365.25, outcome2]/365.25)
  
  # p value is defined by the CC vs. TT comparison, not the dosage
  data$pred2.mod = ifelse(data[,pred2] == 0, 0, ifelse(data[,pred2] == 2, 1, NA))
  p = cmprsk::cuminc(ftime = data[, outcome2],
                     fstatus = data[, status2],
                     group = data[, "pred2.mod"])
  pcomb$p.cuminc2.mod = ifelse(p$Tests[1,2] < 0.0001, "<0.0001", signif(p$Tests[1,2], 2))
  
  #
  
  return(select(pcomb, var, 
                inc.10y.0.1, inc.10y.1.1, inc.10y.2.1, inc.10y.1.diff, p.cuminc1.mod, events.1, fu.1,
                inc.10y.0.2, inc.10y.1.2, inc.10y.2.2, inc.10y.2.diff, p.cuminc2.mod, events.2, fu.2))
  
}


#### load files ####
data0 = as.data.frame(read_csv("file1.csv")) # complete cohort
data = as.data.frame(read_csv("file2.csv")) # prospective cohort

#### descriptive text ####
# In the complete cohort, MACE occurred in X participants at any time
count(data0, status.cv.any == 1) #52024

# and in X by age 70, 
count(data0, status.cv.any == 1 & age.cv.any <= 70) #41765

# while LRE occurred in X and X at any time and by age 70, respectively.
count(data0, status.c == 1) # 1413
count(data0, status.c == 1 & age.cd <= 70) # 1003

# In the prospective cohort, X and X participants developed MACE and LRE, respectively. 
count(data, status.cv.any == 1) #28347
count(data, status.c == 1) #850

data$fu = as.numeric(data$date.last - data$date.enroll)
# Median follow-up in the prospective cohort was X months and 
median(data$fu)/30.5

# total follow-up was X person-years.
sum(data$fu)/365.25


#### Table 1: descriptive ####
d.cat = as.data.frame(select(data0,
                             female,
                             asian, black, white, other,
                             lean, ow, o1, o2, o3,
                             dm.prev, hl.prev, htn.prev,
                             FIB4.1, FIB4.2, FIB4.3,
                             APRI.1, APRI.2, APRI.3))
d.num = as.data.frame(select(data0,
                             agei1,
                             AST1, ALT1, PLT1, ALB1, TG, HDL, LDL, HBA1C))

d.cat1 = as.data.frame(select(dplyr::filter(data, tm6sf2 == 0),
                              female,
                              asian, black, white, other,
                              lean, ow, o1, o2, o3,
                              dm.prev, hl.prev, htn.prev,
                              FIB4.1, FIB4.2, FIB4.3,
                              APRI.1, APRI.2, APRI.3))
d.num1 = as.data.frame(select(dplyr::filter(data, tm6sf2 == 0),
                              agei1,
                              AST1, ALT1, PLT1, ALB1, TG, HDL, LDL, HBA1C))

d.cat2 = as.data.frame(select(dplyr::filter(data, tm6sf2 == 1),
                              female,
                              asian, black, white, other,
                              lean, ow, o1, o2, o3,
                              dm.prev, hl.prev, htn.prev,
                              FIB4.1, FIB4.2, FIB4.3,
                              APRI.1, APRI.2, APRI.3))
d.num2 = as.data.frame(select(dplyr::filter(data, tm6sf2 == 1),
                              agei1,
                              AST1, ALT1, PLT1, ALB1, TG, HDL, LDL, HBA1C))

d.cat3 = as.data.frame(select(dplyr::filter(data, tm6sf2 == 2),
                              female,
                              asian, black, white, other,
                              lean, ow, o1, o2, o3,
                              dm.prev, hl.prev, htn.prev,
                              FIB4.1, FIB4.2, FIB4.3,
                              APRI.1, APRI.2, APRI.3))
d.num3 = as.data.frame(select(dplyr::filter(data, tm6sf2 == 2),
                              agei1,
                              AST1, ALT1, PLT1, ALB1, TG, HDL, LDL, HBA1C))


###### numerical variables ####
table1temp <- names(d.num1)
table1temp <- as.data.frame(table1temp)
table1temp[,1] <- as.character(table1temp[,1])

table1temp$mean1 = sapply(names(d.num1), function(x) mean(d.num1[,x], na.rm = TRUE))        # means
table1temp <- table1temp %>% mutate(se1 = sapply(d.num1, function(...) { sd(..., na.rm = TRUE)/sqrt(length(...)) }),
                                    sd1 = sapply(d.num1, function(...) { sd(..., na.rm = TRUE)}))
table1temp$mean2 = sapply(names(d.num1), function(x) mean(d.num2[,x], na.rm = TRUE))        # means
table1temp <- table1temp %>% mutate(se2 = sapply(d.num2, function(...) { sd(..., na.rm = TRUE)/sqrt(length(...)) }),
                                    sd2 = sapply(d.num2, function(...) { sd(..., na.rm = TRUE)}))
table1temp$mean3 = sapply(names(d.num1), function(x) mean(d.num3[,x], na.rm = TRUE))        # means
table1temp <- table1temp %>% mutate(se3 = sapply(d.num3, function(...) { sd(..., na.rm = TRUE)/sqrt(length(...)) }),
                                    sd3 = sapply(d.num3, function(...) { sd(..., na.rm = TRUE)}))
table1temp$n = sapply(table1temp$table1temp, function(x) sum(!is.na(d.num[,x])))

table1temp$p = sapply(table1temp$table1temp, function(x) {
  summary(aov(as.formula(paste0(x, "~ as.factor(tm6sf2)")), data = data))[[1]][["Pr(>F)"]][1]
})


# formatting it so that it looks nicer
table1 <- NULL
table1$title <- table1temp$table1temp
table1 <- as.data.frame(table1)

# this code below rounds to the number of significant figures that you want
# sprintf("%.1f"... will round to 1 digit (0.1). If you want 2 digits, use: sprintf("%.2f"...
table1$tm6sf2.0 <- paste0(sprintf("%.1f", table1temp$mean1), " (", sprintf("%.1f", table1temp$sd1), ")")
table1$tm6sf2.1 <- paste0(sprintf("%.1f", table1temp$mean2), " (", sprintf("%.1f", table1temp$sd2), ")")
table1$tm6sf2.2 <- paste0(sprintf("%.1f", table1temp$mean3), " (", sprintf("%.1f", table1temp$sd3), ")")

table1$n <- table1temp$n

# formatting P value
table1$p <- ifelse(table1temp$p < 0.0001, "<0.0001", sprintf("%.4f", table1temp$p))

###### categorical variables ####
table1atemp <- names(d.cat)
table1atemp <- as.data.frame(table1atemp)
table1atemp[,1] <- as.character(table1atemp[,1])

table1atemp <- table1atemp %>% 
  mutate(
    
    prop1 = sapply(names(d.cat1), function(x) {
      sum(d.cat1[,x] == 1, na.rm = T) / 
        (sum(d.cat1[,x] == 1, na.rm = T) + sum(d.cat1[,x] == 0, na.rm = T))}),
    
    # n1.yes and n1.no are the NUMBER of people in each category.
    # this is necessary to run a chi-squared test.
    n1.yes = sapply(names(d.cat1), function(x) {
      sum(d.cat1[,x] == 1, na.rm = T)}),
    n1.no = sapply(names(d.cat1), function(x) {
      sum(d.cat1[,x] == 0, na.rm = T)}),
    
    # same for categories 2 and 3
    prop2 = sapply(names(d.cat2), function(x) {
      sum(d.cat2[,x] == 1, na.rm = T) / 
        (sum(d.cat2[,x] == 1, na.rm = T) + sum(d.cat2[,x] == 0, na.rm = T))}),
    n2.yes = sapply(names(d.cat2), function(x) {
      sum(d.cat2[,x] == 1, na.rm = T)}),
    n2.no = sapply(names(d.cat2), function(x) {
      sum(d.cat2[,x] == 0, na.rm = T)}),
    
    prop3 = sapply(names(d.cat3), function(x) {
      sum(d.cat3[,x] == 1, na.rm = T) / 
        (sum(d.cat3[,x] == 1, na.rm = T) + sum(d.cat3[,x] == 0, na.rm = T))}),
    n3.yes = sapply(names(d.cat3), function(x) {
      sum(d.cat3[,x] == 1, na.rm = T)}),
    n3.no = sapply(names(d.cat3), function(x) {
      sum(d.cat3[,x] == 0, na.rm = T)})
  )

table1atemp$n = sapply(table1atemp$table1atemp, function(x) sum(!is.na(d.cat[,x])))

table1atemp$p = sapply(table1atemp$table1atemp, function(x) {
  chisq.test(matrix(c(table1atemp[table1atemp$table1atemp == x, 3],
                      table1atemp[table1atemp$table1atemp == x, 4],
                      table1atemp[table1atemp$table1atemp == x, 6],
                      table1atemp[table1atemp$table1atemp == x, 7],
                      table1atemp[table1atemp$table1atemp == x, 9],
                      table1atemp[table1atemp$table1atemp == x, 10]), 
                    nrow = 2, ncol = 3, byrow = T))$p.value
})

# again just cleaning up the variable, similar to as above
table1a <- NULL
table1a$title <- table1atemp$table1atemp
table1a <- as.data.frame(table1a)

# convert proportions into percentages
table1a$tm6sf2.0 <- paste(as.character(sprintf("%.1f", table1atemp$prop1 * 100)), "%", sep = "")
table1a$tm6sf2.1 <- paste(as.character(sprintf("%.1f", table1atemp$prop2 * 100)), "%", sep = "")
table1a$tm6sf2.2 <- paste(as.character(sprintf("%.1f", table1atemp$prop3 * 100)), "%", sep = "")

table1a$n <- table1atemp$n
table1a$p <- ifelse(table1atemp$p < 0.0001, "<0.0001", sprintf("%.4f", table1atemp$p))

# combine the "cleaned" 2 tables and you're done
table.final = rbind(table1, table1a)

write_csv(table.final, "TABLE1.csv")




#### Table 2: Prevalence of liver-related events and major adverse cardiovascular events by age 70 years based on TM6SF2 genotype ####

hosp.age = function(data, lab = "") {
  
  x1 = count(data %>% dplyr::filter(!is.na(age.lre), !is.na(tm6sf2)), age.lre > 70, round(tm6sf2), status.lre == 1)
  x2 = count(data %>% dplyr::filter(!is.na(age.mace), !is.na(tm6sf2)), age.mace > 70, round(tm6sf2), status.mace == 1)
  
  x1.0.pos = ifelse(length(x1[x1[,1] == F & x1[,3] == T & x1[,2] == 0, "n"]) == 1,
                    x1[x1[,1] == F & x1[,3] == T & x1[,2] == 0, "n"], 0)
  x1.1.pos = ifelse(length(x1[x1[,1] == F & x1[,3] == T & x1[,2] == 1, "n"]) == 1,
                    x1[x1[,1] == F & x1[,3] == T & x1[,2] == 1, "n"], 0)
  x1.2.pos = ifelse(length(x1[x1[,1] == F & x1[,3] == T & x1[,2] == 2, "n"]) == 1,
                    x1[x1[,1] == F & x1[,3] == T & x1[,2] == 2, "n"], 0)
  
  x1.0.neg = ifelse(length(x1[x1[,1] == T & x1[,3] == F & x1[,2] == 0, "n"]) == 1,
                    x1[x1[,1] == T & x1[,3] == F & x1[,2] == 0, "n"], 0)
  x1.1.neg = ifelse(length(x1[x1[,1] == T & x1[,3] == F & x1[,2] == 1, "n"]) == 1,
                    x1[x1[,1] == T & x1[,3] == F & x1[,2] == 1, "n"], 0)
  x1.2.neg = ifelse(length(x1[x1[,1] == T & x1[,3] == F & x1[,2] == 2, "n"]) == 1,
                    x1[x1[,1] == T & x1[,3] == F & x1[,2] == 2, "n"], 0)
  
  x2.0.pos = ifelse(length(x2[x2[,1] == F & x2[,3] == T & x2[,2] == 0, "n"]) == 1,
                    x2[x2[,1] == F & x2[,3] == T & x2[,2] == 0, "n"], 0)
  x2.1.pos = ifelse(length(x2[x2[,1] == F & x2[,3] == T & x2[,2] == 1, "n"]) == 1,
                    x2[x2[,1] == F & x2[,3] == T & x2[,2] == 1, "n"], 0)
  x2.2.pos = ifelse(length(x2[x2[,1] == F & x2[,3] == T & x2[,2] == 2, "n"]) == 1,
                    x2[x2[,1] == F & x2[,3] == T & x2[,2] == 2, "n"], 0)
  
  x2.0.neg = ifelse(length(x2[x2[,1] == T & x2[,3] == F & x2[,2] == 0, "n"]) == 1,
                    x2[x2[,1] == T & x2[,3] == F & x2[,2] == 0, "n"], 0)
  x2.1.neg = ifelse(length(x2[x2[,1] == T & x2[,3] == F & x2[,2] == 1, "n"]) == 1,
                    x2[x2[,1] == T & x2[,3] == F & x2[,2] == 1, "n"], 0)
  x2.2.neg = ifelse(length(x2[x2[,1] == T & x2[,3] == F & x2[,2] == 2, "n"]) == 1,
                    x2[x2[,1] == T & x2[,3] == F & x2[,2] == 2, "n"], 0)
  
  table = NULL
  table$lab = lab
  table$age = i
  table = as.data.frame(table)
  
  table$lre.proc.0 = 100*x1.0.pos / (x1.0.pos + x1.0.neg)
  table$lre.proc.1 = 100*x1.1.pos / (x1.1.pos + x1.1.neg)
  table$lre.proc.2 = 100*x1.2.pos / (x1.2.pos + x1.2.neg)
  table$lre.p = chisq.test(matrix(c(x1.0.pos, x1.0.neg, x1.1.pos, x1.1.neg, x1.2.pos, x1.2.neg), nrow = 2))$p.value
  
  table$mace.proc.0 = 100*x2.0.pos / (x2.0.pos + x2.0.neg)
  table$mace.proc.1 = 100*x2.1.pos / (x2.1.pos + x2.1.neg)
  table$mace.proc.2 = 100*x2.2.pos / (x2.2.pos + x2.2.neg)
  table$mace.p = chisq.test(matrix(c(x2.0.pos, x2.0.neg, x2.1.pos, x2.1.neg, x2.2.pos, x2.2.neg), nrow = 2))$p.value
  
  table = mutate(table,
                 
                 lre.0 = paste0(prettyNum(x1.0.pos, big.mark = ","), 
                                 "/",
                                 prettyNum(x1.0.pos + x1.0.neg, big.mark = ","),
                                 " (",
                                 formatC(lre.proc.0, digits = 1, format = "f"), "%)"),
                 lre.1 = paste0(prettyNum(x1.1.pos, big.mark = ","), 
                                 "/",
                                 prettyNum(x1.1.pos + x1.1.neg, big.mark = ","),
                                 " (",
                                 formatC(lre.proc.1, digits = 1, format = "f"), "%)"),
                 lre.2 = paste0(prettyNum(x1.2.pos, big.mark = ","), 
                                 "/",
                                 prettyNum(x1.2.pos + x1.2.neg, big.mark = ","),
                                 " (",
                                 formatC(lre.proc.2, digits = 1, format = "f"), "%)"),
                 
                 mace.0 = paste0(prettyNum(x2.0.pos, big.mark = ","), 
                                "/",
                                prettyNum(x2.0.pos + x2.0.neg, big.mark = ","),
                                " (",
                                formatC(mace.proc.0, digits = 1, format = "f"), "%)"),
                 mace.1 = paste0(prettyNum(x2.1.pos, big.mark = ","),
                                "/",
                                prettyNum(x2.1.pos + x2.1.neg, big.mark = ","),
                                " (",
                                formatC(mace.proc.1, digits = 1, format = "f"), "%)"),
                 mace.2 = paste0(prettyNum(x2.2.pos, big.mark = ","), 
                                "/",
                                prettyNum(x2.2.pos + x2.2.neg, big.mark = ","),
                                " (",
                                formatC(mace.proc.2, digits = 1, format = "f"), "%)"))
  
  return(table %>% select(lab, lre.0, lre.1, lre.2, lre.p, mace.0, mace.1, mace.2, mace.p))
  
}


table2 = hosp.age(data0, lab = "All")
table2 = rbind(table2,
               hosp.age(dplyr::filter(data0, sex == 0), lab = "female"))
table2 = rbind(table2,
               hosp.age(dplyr::filter(data0, sex == 1), lab = "male"))
table2 = rbind(table2, 
               hosp.age(dplyr::filter(data0, dm.prev == 1), lab = "Diabetes"))
table2 = rbind(table2,
               hosp.age(dplyr::filter(data0, FIB4 < 1.3), lab = "FIB4 < 1.3"))
table2 = rbind(table2, 
               hosp.age(dplyr::filter(data0, FIB4 >= 1.3, FIB4 <= 2.67), lab = "FIB4 1.3-2.67"))
table2 = rbind(table2, 
               hosp.age(dplyr::filter(data0, FIB4 > 2.67), lab = "FIB4 > 2.67"))
write_csv(table2, "TABLE2.csv")


#### Table 3: Impact of TM6SF2 genotype on incidence of liver-related events and major adverse cardiovascular events ####
#### and Supp Table 6: E values for associations of TM6SF2 genotype and liver-related events or major adverse cardiovascular events ####

# unadjusted
file = "T3 unadjusted.csv"

table1 = NULL
table1 = rbind(table1,
               hr.table(t = data,
                        rowname = "All"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 0), 
                        rowname = "females"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 1), 
                        rowname = "males"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, dm.prev == T), 
                        rowname = "DM"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 < 1.3),
                        rowname = "FIB4 low"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 >= 1.3, FIB4 <= 2.67),
                        rowname = "FIB4 int"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 > 2.67),
                        rowname = "FIB4 high"))
write_csv(table1, file)

# adjusted
cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1"
file = "T3 adjusted.csv"

# all
table1 = NULL
table1 = rbind(table1,
               hr.table(t = data,
                        cov = cov,
                        rowname = "All"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 1), 
                        cov =  "age+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",
                        rowname = "males"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 0), 
                        cov =  "age+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",
                        rowname = "females"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, dm.prev == T), 
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + htn.prev + alc.high + smoking + statin1",
                        rowname = "DM"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 < 1.3),
                        cov = cov,
                        rowname = "FIB4 low"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 >= 1.3, FIB4 <= 2.67),
                        cov = cov,
                        rowname = "FIB4 int"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 > 2.67),
                        cov = cov,
                        rowname = "FIB4 high"))
write_csv(table1, file)



#### Table 4: 10 year cumulative incidence of liver-related events and major adverse cardiovascular events based on TM6SF2 genotype ####
tab = cuminc10y(data = data, 
                rowname = "All")
tab = rbind(tab,
            cuminc10y(data = dplyr::filter(data, sex == 0), 
                      rowname = "Females"))
tab = rbind(tab,
            cuminc10y(data = dplyr::filter(data, sex == 1), 
                      rowname = "Males"))
tab = rbind(tab,
            cuminc10y(data = dplyr::filter(data, DM == 1), 
                      rowname = "Diabetes"))
tab = rbind(tab,
            cuminc10y(data = dplyr::filter(data, FIB4 < 1.3), 
                      rowname = "FIB4 low"))
tab = rbind(tab,
            cuminc10y(data = dplyr::filter(data, FIB4 >= 1.3, FIB4 <= 2.67), 
                      rowname = "FIB4 int"))
tab = rbind(tab,
            cuminc10y(data = dplyr::filter(data, FIB4 > 2.67), 
                      rowname = "FIB4 high"))
write_csv(tab, "T4.csv")



#### Table 5: Unadjusted associations between TM6SF2 genotype and liver magnetic resonance imaging properties ####
summ.reg = function(regression, type = "lin") {
  a = summary(regression)$coef
  
  b=NULL
  b$var = rownames(a)
  b = as.data.frame(b)
  b$beta = a[,1]
  b$se = a[,2]
  #b$p = format(a[, ncol(a)], scientific = T, digits = 2)
  b$p = a[, ncol(a)]
  
  if (type == "log") {
    b = mutate(b, summ = paste0(sprintf("%.2f", exp(beta)), " (",
                                sprintf("%.2f", exp(beta - 1.96*se)), "-",
                                sprintf("%.2f", exp(beta + 1.96*se)), ")"))
  } else if (type == "lin") {
    b = mutate(b, summ = paste0(sprintf("%.2f", beta), " (",
                                sprintf("%.2f", beta - 1.96 * se), " to ",
                                sprintf("%.2f", beta + 1.96 * se), ")"))
  }
  return(select(b, var, beta, se, summ, p))
}


uv.analysis2 = function(data1, cov, name, date) {
  
  outcomes = c("PDFF", "cT1", "PDFF.h", "cT1.h")
  j = "sex"
  cov.all = base::strsplit(cov, "[ ]*\\+[ ]*")[[1]]
  
  for (outcome in outcomes) {
    for (k in 0:1) {
      if (outcome %in% c("PDFF", "cT1")) {
        a = summ.reg(glm(as.formula(paste0(outcome," ~ ", cov, "tm6sf2")), data = data1[data1[,j] == k,]), type = "lin")
        write_csv(a, paste0(name, "-", outcome, "-", j, "-", k, "-", date, ".csv"))
      } else if (outcome %in% c("PDFF.h", "cT1.h")) {
        a = summ.reg(glm(as.formula(paste0(outcome," ~ ", cov, "tm6sf2")), data = data1[data1[,j] == k,]), type = "log")
        write_csv(a, paste0(name, "-", outcome, "-", j, "-", k, "-", date, ".csv"))
      }
    }
  }
  
  table = NULL
  table$name = outcomes
  table = as.data.frame(table)
  for (outcome in outcomes) {
    for (k in 0:1) {
      a = as.data.frame(read_csv(show_col_types = FALSE, paste0(name, "-", outcome, "-", j, "-", k, "-", date, ".csv")))
      table[table$name == outcome, paste0("summ.", k)] = a[a$var == "tm6sf2", "summ"]
      table[table$name == outcome, paste0("p.", k)] = a[a$var == "tm6sf2", "p"]
    }
    x = hetero.uv1(outcome = outcome, name = name, j = j, 
                   k1 = 0, k2 = 1, date = date)
    table[table$name == outcome, "phet"] = x[1, "phet"]
    table[table$name == outcome, "n"] = dim(data1[complete.cases(data1[, c(cov.all, "tm6sf2", outcome)]), ])[1]
  }
  
  write_csv(table, paste0(name, "-", j, "-", date, ".csv"))
  return(table)
}


hetero.uv1 = function(outcome = "PDFF", name, j, k1, k2, date) {
  
  t2 = NULL
  t2$diet = outcome
  t2 = as.data.frame(t2)
  
  data1 = read_csv(show_col_types = FALSE, paste0(name, "-", outcome, "-", j, "-", k1, "-", date, ".csv"))
  data2 = read_csv(show_col_types = FALSE, paste0(name, "-", outcome, "-", j, "-", k2, "-", date, ".csv"))
  
  t2[t2$diet == outcome, "beta1"] = data1[data1$var == "tm6sf2", "beta"]
  t2[t2$diet == outcome, "se1"] = data1[data1$var == "tm6sf2", "se"]
  t2[t2$diet == outcome, "beta2"] = data2[data2$var == "tm6sf2", "beta"]
  t2[t2$diet == outcome, "se2"] = data2[data2$var == "tm6sf2", "se"]
  
  t2 = mutate(t2,
              beta.pooled = 
                (beta1/(se1^2) + beta2/(se2^2))/
                (1/se1^2 + 1/se2^2),
              q = 
                (beta1 - beta.pooled)^2*(1/se1^2) +
                (beta2 - beta.pooled)^2*(1/se2^2),
              phet = pchisq(q, df = 1, lower.tail = F))
  
  return(select(t2, diet, phet))
  
}

date = "2024-08-30"
cov=""
a = uv.analysis2(data1 = data0,
                 cov = cov, name = "unadj-all", date = date)
a$cat = "unadj-all"

b = uv.analysis2(data1 = data0 %>% filter(FIB4 < 1.3),
                 cov = cov, name = "unadj-FIB4.low", date = date)
b$cat = "unadj-FIB4.low"

c = uv.analysis2(data1 = data0 %>% filter(FIB4 >= 1.3, FIB4 <= 2.67),
                 cov = cov, name = "unadj-FIB4.int", date = date)
c$cat = "unadj-FIB4.int"

d = uv.analysis2(data1 = data0 %>% filter(FIB4 > 2.67),
                 cov = cov, name = "unadj-FIB4.high", date = date)
d$cat = "unadj-FIB4.high"

e = uv.analysis2(data1 = data0 %>% filter(DM == 1),
                 cov = cov, name = "unadj-DM", date = date)
e$cat = "unadj-DM"

table = rbind(a, b, c, d, e)
write_csv(table[, c(8, 1:7)], "T5.csv")


#### Supp Table 3: Liver-related death in the overall cohort by age 70 ####

data0$death.age = (data0$death.date - data0$dob)/365.25

death.age = function(data, lab = "") {
  i=70
  x1 = count(data %>% dplyr::filter(!is.na(age.lre), !is.na(tm6sf2)),
             death.age < i & deaths.liv == T, round(tm6sf2))
  x2 = count(data %>% dplyr::filter(!is.na(age.mace), !is.na(tm6sf2)), 
             death.age < i & deaths.mace == T, round(tm6sf2))
  
  x1.0.pos = ifelse(length(x1[x1[,1] == T & x1[,2] == 0, "n"]) == 1,
                    x1[x1[,1] == T & x1[,2] == 0, "n"], 0)
  x1.1.pos = ifelse(length(x1[x1[,1] == T & x1[,2] == 1, "n"]) == 1,
                    x1[x1[,1] == T & x1[,2] == 1, "n"], 0)
  x1.2.pos = ifelse(length(x1[x1[,1] == T & x1[,2] == 2, "n"]) == 1,
                    x1[x1[,1] == T & x1[,2] == 2, "n"], 0)
  
  x1.0.neg = ifelse(length(x1[x1[,1] == F & x1[,2] == 0, "n"]) == 1,
                    x1[x1[,1] == F & x1[,2] == 0, "n"], 0)
  x1.1.neg = ifelse(length(x1[x1[,1] == F & x1[,2] == 1, "n"]) == 1,
                    x1[x1[,1] == F & x1[,2] == 1, "n"], 0)
  x1.2.neg = ifelse(length(x1[x1[,1] == F &  x1[,2] == 2, "n"]) == 1,
                    x1[x1[,1] == F & x1[,2] == 2, "n"], 0)
  
  
  x2.0.pos = ifelse(length(x2[x2[,1] == T & x2[,2] == 0, "n"]) == 1,
                    x2[x2[,1] == T & x2[,2] == 0, "n"], 0)
  x2.1.pos = ifelse(length(x2[x2[,1] == T & x2[,2] == 1, "n"]) == 1,
                    x2[x2[,1] == T & x2[,2] == 1, "n"], 0)
  x2.2.pos = ifelse(length(x2[x2[,1] == T & x2[,2] == 2, "n"]) == 1,
                    x2[x2[,1] == T & x2[,2] == 2, "n"], 0)
  
  x2.0.neg = ifelse(length(x2[x2[,1] == F & x2[,2] == 0, "n"]) == 1,
                    x2[x2[,1] == F & x2[,2] == 0, "n"], 0)
  x2.1.neg = ifelse(length(x2[x2[,1] == F & x2[,2] == 1, "n"]) == 1,
                    x2[x2[,1] == F & x2[,2] == 1, "n"], 0)
  x2.2.neg = ifelse(length(x2[x2[,1] == F &  x2[,2] == 2, "n"]) == 1,
                    x2[x2[,1] == F & x2[,2] == 2, "n"], 0)
  
  table = NULL
  table$lab = lab
  table$age = i
  table = as.data.frame(table)
  table$lre.proc.0 = 100*x1.0.pos / (x1.0.pos + x1.0.neg)
  table$lre.proc.1 = 100*x1.1.pos / (x1.1.pos + x1.1.neg)
  table$lre.proc.2 = 100*x1.2.pos / (x1.2.pos + x1.2.neg)
  table$lre.p = chisq.test(matrix(c(x1.0.pos, x1.0.neg, x1.1.pos, x1.1.neg, x1.2.pos, x1.2.neg), nrow = 2))$p.value
  table$mace.proc.0 = 100*x2.0.pos / (x2.0.pos + x2.0.neg)
  table$mace.proc.1 = 100*x2.1.pos / (x2.1.pos + x2.1.neg)
  table$mace.proc.2 = 100*x2.2.pos / (x2.2.pos + x2.2.neg)
  table$mace.p = chisq.test(matrix(c(x2.0.pos, x2.0.neg, x2.1.pos, x2.1.neg, x2.2.pos, x2.2.neg), nrow = 2))$p.value
  
  table = mutate(table,
                 
                 lre.0 = paste0(prettyNum(x1.0.pos, big.mark = ","), 
                                "/",
                                prettyNum(x1.0.pos + x1.0.neg, big.mark = ","),
                                " (",
                                formatC(lre.proc.0, digits = 1, format = "f"), "%)"),
                 lre.1 = paste0(prettyNum(x1.1.pos, big.mark = ","), 
                                "/",
                                prettyNum(x1.1.pos + x1.1.neg, big.mark = ","),
                                " (",
                                formatC(lre.proc.1, digits = 1, format = "f"), "%)"),
                 lre.2 = paste0(prettyNum(x1.2.pos, big.mark = ","), 
                                "/",
                                prettyNum(x1.2.pos + x1.2.neg, big.mark = ","),
                                " (",
                                formatC(lre.proc.2, digits = 1, format = "f"), "%)"),
                 
                 mace.0 = paste0(prettyNum(x2.0.pos, big.mark = ","), 
                                 "/",
                                 prettyNum(x2.0.pos + x2.0.neg, big.mark = ","),
                                 " (",
                                 formatC(mace.proc.0, digits = 1, format = "f"), "%)"),
                 mace.1 = paste0(prettyNum(x2.1.pos, big.mark = ","),
                                 "/",
                                 prettyNum(x2.1.pos + x2.1.neg, big.mark = ","),
                                 " (",
                                 formatC(mace.proc.1, digits = 1, format = "f"), "%)"),
                 mace.2 = paste0(prettyNum(x2.2.pos, big.mark = ","), 
                                 "/",
                                 prettyNum(x2.2.pos + x2.2.neg, big.mark = ","),
                                 " (",
                                 formatC(mace.proc.2, digits = 1, format = "f"), "%)"))
  
  return(table %>% select(lab, lre.0, lre.1, lre.2, lre.p, mace.0, mace.1, mace.2, mace.p))
  
}


ST3 = death.age(data0, lab = "All")
ST3 = rbind(ST3,
            death.age(dplyr::filter(data0, sex == 0), lab = "female"))
ST3 = rbind(ST3,
            death.age(dplyr::filter(data0, sex == 1), lab = "male"))
ST3 = rbind(ST3, 
            death.age(dplyr::filter(data0, DM == 1), lab = "Diabetes"))
ST3 = rbind(ST3,
            death.age(dplyr::filter(data0, FIB4 < 1.3), lab = "FIB4 < 1.3"))
ST3 = rbind(ST3, 
            death.age(dplyr::filter(data0, FIB4 >= 1.3, FIB4 <= 2.67), lab = "FIB4 1.3-2.67"))
ST3 = rbind(ST3, 
            death.age(dplyr::filter(data0, FIB4 > 2.67), lab = "FIB4 > 2.67"))
write_csv(table2, "SUPP TABLE 3.csv")


#### Supp Table 4: Subgroup analyses of associations between TM6SF2 genotype and major adverse cardiovascular events or liver related events ####

# unadjusted
file = "ST4 unadjusted.csv"

table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, round(rs738409) == 0),
                        rowname = "pnpla3 wt"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, round(rs738409) > 0),
                        rowname = "pnpla3 abnl"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, anc == "EUR"),
                        rowname = "eur"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, anc != "EUR"),
                        rowname = "noneur"))
write_csv(table1, file)

# adjusted
cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1"
file = "ST4 adjusted.csv"
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, round(rs738409) == 0),
                        cov = cov,
                        rowname = "pnpla3 wt"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, round(rs738409) > 0),
                        cov = cov,
                        rowname = "pnpla3 abnl"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, anc == "EUR"),
                        cov = cov,
                        rowname = "eur"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, anc != "EUR"),
                        cov = cov,
                        rowname = "noneur"))
write_csv(table1, file)




#### Supp. Table 5: Impact of TM6SF2 genotype on incidence of liver- and cardiac-related mortality ####
file = "ST5 unadjusted.csv"
table1 = hr.table(t = data,
                  time1 = 'time.death.lre', status1 = 'status.death.lre',
                  time2 = "time.death.mace", status2 = 'status.death.mace',
                  rowname = "All")
write_csv(table1, file)

file = "ST5 adjusted.csv"
cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1"
table1 = hr.table(t = data,
                  time1 = 'time.death.lre', status1 = 'status.death.lre',
                  time2 = "time.death.mace", status2 = 'status.death.mace',
                  cov = cov,
                  rowname = "All")
write_csv(table1, file)


#### Supp Table 6: see code for Table 3 ####

#### Supp Table 7: Associations between TM6SF2 genotype and major adverse cardiovascular events or liver related events, stratified by sex and Fibrosis-4 score ####

# unadjusted

file = "ST7 unadjusted.csv"

table1 = NULL
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 < 1.3, sex == 1),
                        rowname = "FIB4 low male"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 >= 1.3, FIB4 <= 2.67, sex == 1),
                        rowname = "FIB4 int male"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 > 2.67, sex == 1),
                        rowname = "FIB4 high male"))
write_csv(table1, file)

table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 < 1.3, sex == 0),
                        rowname = "FIB4 low female"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 >= 1.3, FIB4 <= 2.67, sex == 0),
                        rowname = "FIB4 int female"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, FIB4 > 2.67, sex == 0),
                        rowname = "FIB4 high female"))
write_csv(table1, file)

# adjusted
file = "ST7 adjusted.csv"
cov = "age+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1"

table1 = NULL
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 1, FIB4 < 1.3),
                        cov = cov,
                        rowname = "FIB4 low male"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 1, FIB4 >= 1.3, FIB4 <= 2.67),
                        cov = cov,
                        rowname = "FIB4 int male"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 1, FIB4 > 2.67),
                        cov = cov,
                        rowname = "FIB4 high male"))
write_csv(table1, file)

table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 0, FIB4 < 1.3),
                        cov = cov,
                        rowname = "FIB4 low female"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 0, FIB4 >= 1.3, FIB4 <= 2.67),
                        cov = cov,
                        rowname = "FIB4 int female"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = dplyr::filter(data, sex == 0, FIB4 > 2.67),
                        cov = cov,
                        rowname = "FIB4 high female"))
write_csv(table1, file)

#### Supp. Table 8: Effects of increasing landmark periods times on associations between TM6SF2 genotype and liver-related events and major adverse cardiovascular events ####
# status.mace.1, status.mace.2, status.mace.5, status.lre.1, status.lre.2, status.lre.5 are set to NA if time.mace/time.lre are under 1, 2, and 5 years, respectively

# unadjusted
file = "ST8 unadjusted.csv"
table1 = hr.table(t = data,
                  status1 = 'status.mace.1', time1 = "time.mace",
                  status2 = "status.lre.1", time2 = "time.lre",
                  rowname = "1yr")
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data,
                        status1 = 'status.mace.2', time1 = "time.mace",
                        status2 = "status.lre.2", time2 = "time.lre",
                        rowname = "2yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data,
                        status1 = 'status.mace.5', time1 = "time.mace",
                        status2 = "status.lre.5", time2 = "time.lre",
                        rowname = "5yr"))
write_csv(table1, file)

table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 < 1.3),
                        status1 = 'status.mace.1', time1 = "time.mace",
                        status2 = "status.lre.1", time2 = "time.lre",
                        rowname = "FIB4 low - 1yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 < 1.3),
                        status1 = 'status.mace.2', time1 = "time.mace",
                        status2 = "status.lre.2", time2 = "time.lre",
                        rowname = "FIB4 low - 2yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 < 1.3),
                        status1 = 'status.mace.5', time1 = "time.mace",
                        status2 = "status.lre.5", time2 = "time.lre",
                        rowname = "FIB4 low - 5yr"))
write_csv(table1, file)

table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 >= 1.3, FIB4 <= 2.67),
                        status1 = 'status.mace.1', time1 = "time.mace",
                        status2 = "status.lre.1", time2 = "time.lre",
                        rowname = "FIB4 int - 1yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 >= 1.3, FIB4 <= 2.67),
                        status1 = 'status.mace.2', time1 = "time.mace",
                        status2 = "status.lre.2", time2 = "time.lre",
                        rowname = "FIB4 int - 2yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 >= 1.3, FIB4 <= 2.67),
                        status1 = 'status.mace.5', time1 = "time.mace",
                        status2 = "status.lre.5", time2 = "time.lre",
                        rowname = "FIB4 int - 5yr"))
write_csv(table1, file)

table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 > 2.67),
                        status1 = 'status.mace.1', time1 = "time.mace",
                        status2 = "status.lre.1", time2 = "time.lre",
                        rowname = "FIB4 high - 1yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 > 2.67),
                        status1 = 'status.mace.2', time1 = "time.mace",
                        status2 = "status.lre.2", time2 = "time.lre",
                        rowname = "FIB4 high - 2yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 >2.67),
                        status1 = 'status.mace.5', time1 = "time.mace",
                        status2 = "status.lre.5", time2 = "time.lre",
                        rowname = "FIB4 high - 5yr"))
write_csv(table1, file)

# adjusted

file = "ST8 adjusted.csv"

table1 = hr.table(t = data,
                  cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",
                  status1 = 'status.mace.1', time1 = "time.mace",
                  status2 = "status.lre.1", time2 = "time.lre",
                  rowname = "1yr")
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data,
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",
                        status1 = 'status.mace.2', time1 = "time.mace",
                        status2 = "status.lre.2", time2 = "time.lre",
                        rowname = "2yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data,
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",
                        status1 = 'status.mace.5', time1 = "time.mace",
                        status2 = "status.lre.5", time2 = "time.lre",
                        rowname = "5yr"))
write_csv(table1, file)

table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 < 1.3),
                        status1 = 'status.mace.1', time1 = "time.mace",
                        status2 = "status.lre.1", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",  
                        rowname = "FIB4 low - 1yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 < 1.3),
                        status1 = 'status.mace.2', time1 = "time.mace",
                        status2 = "status.lre.2", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",  
                        rowname = "FIB4 low - 2yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 < 1.3),
                        status1 = 'status.mace.5', time1 = "time.mace",
                        status2 = "status.lre.5", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",  
                        rowname = "FIB4 low - 5yr"))
write_csv(table1, file)

table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 >= 1.3, FIB4 <= 2.67),
                        status1 = 'status.mace.1', time1 = "time.mace",
                        status2 = "status.lre.1", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",  
                        rowname = "FIB4 int - 1yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 >= 1.3, FIB4 <= 2.67),
                        status1 = 'status.mace.2', time1 = "time.mace",
                        status2 = "status.lre.2", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",  
                        rowname = "FIB4 int - 2yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 >= 1.3, FIB4 <= 2.67),
                        status1 = 'status.mace.5', time1 = "time.mace",
                        status2 = "status.lre.5", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",  
                        rowname = "FIB4 int - 5yr"))
write_csv(table1, file)

table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 > 2.67),
                        status1 = 'status.mace.1', time1 = "time.mace",
                        status2 = "status.lre.1", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1", 
                        rowname = "FIB4 high - 1yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 > 2.67),
                        status1 = 'status.mace.2', time1 = "time.mace",
                        status2 = "status.lre.2", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1", 
                        rowname = "FIB4 high - 2yr"))
write_csv(table1, file)
table1 = rbind(table1,
               hr.table(t = data %>% filter(FIB4 >2.67),
                        status1 = 'status.mace.5', time1 = "time.mace",
                        status2 = "status.lre.5", time2 = "time.lre",
                        cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",  
                        rowname = "FIB4 high - 5yr"))
write_csv(table1, file)

#### Supp. Table 9: Impact of TM6SF2 genotype on incidence of liver-related events based on an additive model and major adverse cardiovascular events based on a recessive model ####
hr.table = function(t, rowname = "", time1 = "time.mace", status1 = "status.mace", time2 = "time.lre", status2 = "status.lre", cov) {
  
  # for MACE outcome (status1), set recessive vs. not
  # for LRE outcome (status2), set additive model - tm6sf2 is an integer 0-2
  if (cov == "") {
    cov.tot1 = "(tm6sf2 == 2)"
    cov.tot2 = "tm6sf2"
  } else {
    cov.tot1 = paste("(tm6sf2 == 2)", cov, sep = "+")
    cov.tot2 = paste("tm6sf2", cov, sep = "+")
  }
  cov.all = base::strsplit(cov.tot2, "[ ]*\\+[ ]*")[[1]]
  
  # recode status variable to ensure consistent with tidycmprsk
  t[ , status1] = as.factor(case_when(t[, status1] == 1 ~ 2,
                                      t[, status1] == 2 ~ 3,
                                      t[, status1] == 3 ~ 1))
  
  t[ , status2] = as.factor(case_when(t[, status2] == 1 ~ 2,
                                      t[, status2] == 2 ~ 3,
                                      t[, status2] == 3 ~ 1))
  
  # convert from days to years
  t[ , time1] = t[, time1]/365.25
  t[ , time2] = t[, time2]/365.25
  
  # run models
  a = crr(as.formula(paste0("Surv(", time1, ",", status1, ") ~ ", cov.tot)), data = t, failcode = 2)
  b = crr(as.formula(paste0("Surv(", time2, ",", status2, ") ~ ", cov.tot)), data = t, failcode = 2)
  
  # extract data from models
  a1 = as.data.frame(a$tidy)
  b1 = as.data.frame(b$tidy)
  
  table = data.frame(
    name = rep(rowname, 2),
    
    HR1 = c(paste0(sprintf(exp(a1[1,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(a1[1,2] - 1.96*a1[1,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(a1[1,2] + 1.96*a1[1,3]), fmt = "%#.2f"), ")"),
            paste0(sprintf(exp(a1[2,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(a1[2,2] - 1.96*a1[2,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(a1[2,2] + 1.96*a1[2,3]), fmt = "%#.2f"), ")")),
    p1 = signif(c(a1[1,5], a1[2,5]), 3),
    
    HR2 = c(paste0(sprintf(exp(b1[1,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(b1[1,2] - 1.96*b1[1,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(b1[1,2] + 1.96*b1[1,3]), fmt = "%#.2f"), ")"),
            paste0(sprintf(exp(b1[2,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(b1[2,2] - 1.96*b1[2,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(b1[2,2] + 1.96*b1[2,3]), fmt = "%#.2f"), ")")),
    p2 = signif(c(b1[1,5], b1[2,5]), 2),
    n = dim(t[rowSums(is.na(t[, cov.all]))==0, ])[1],
    cases.1 = dim(t[which(rowSums(is.na(t[, cov.all]))==0 & t[, status1] == 2), ])[1],
    cases.2 = dim(t[which(rowSums(is.na(t[, cov.all]))==0 & t[, status2] == 2), ])[1])
  
  hr1=exp(a1[2,2])
  ll1=exp(a1[2,2]-1.96*a1[2,3])
  ul1=exp(a1[2,2]+1.96*a1[2,3])
  
  hr2=exp(b1[2,2])
  ll2=exp(b1[2,2]-1.96*b1[2,3])
  ul2=exp(b1[2,2]+1.96*b1[2,3])
  
  e1 = ifelse((hr1 > 1 & ll1 < 1) | # if beta positive but LL negative -> E=1
                (hr1 < 1 & ul1 > 1) | # if beta negative but UL positive -> E=1
                hr1 == 1, 1,  # if HR == 1 -> E=1
              ifelse(hr1 > 1, 
                     ll1 + sqrt(ll1 * (ll1 - 1)),
                     1/ul1 + sqrt(1/ul1 * (1/ul1 - 1))))
  
  e2 = ifelse((hr2 > 1 & ll2 < 1) | # if beta positive but LL negative -> E=1
                (hr2 < 1 & ul2 > 1) | # if beta negative but UL positive -> E=1
                hr2 == 1, 1,  # if HR == 1 -> E=1
              ifelse(hr2 > 1, 
                     ll2 + sqrt(ll2 * (ll2 - 1)),
                     1/ul2 + sqrt(1/ul2 * (1/ul2 - 1))))
  
  table = mutate(table,
                 E1 = e1,
                 E2 = e2)
  
  return(table)
  
}




#### Supp. Table 10: Interactions between TM6SF2 genotype and environmental factors on incidence of liver-related events ####

hr.table.int = function(t, rowname = "", pred, time1 = "time.mace", status1 = "status.mace", time2 = "time.lre", status2 = "status.lre", cov = "") {
  
  x = Sys.time()
  
  t[ , status1] = as.factor(case_when(t[, status1] == 1 ~ 2,
                                      t[, status1] == 2 ~ 3,
                                      t[, status1] == 3 ~ 1))
  
  t[ , status2] = as.factor(case_when(t[, status2] == 1 ~ 2,
                                      t[, status2] == 2 ~ 3,
                                      t[, status2] == 3 ~ 1))
  
  t[ , time1] = t[, time1]/365.25
  t[ , time2] = t[, time2]/365.25
  
  t$t1 = round(t$tm6sf2) == 2
  t$t2 = round(t$tm6sf2)
  
  a = tidycmprsk::crr(stats::as.formula(paste0("Surv(", time1, ",", status1, ") ~ t1*", pred, cov)), data = t, failcode = 2)
  b = tidycmprsk::crr(stats::as.formula(paste0("Surv(", time2, ",", status2, ") ~ t2*", pred, cov)), data = t, failcode = 2)
  
  a1 = as.data.frame(a$tidy)
  b1 = as.data.frame(b$tidy)
  
  table = data.frame(
    name = c("tm6sf2", rowname, paste0("tm6sf2*", rowname)),
    
    HR1 = c(paste0(sprintf(exp(a1[1,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(a1[1,2] - 1.96*a1[1,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(a1[1,2] + 1.96*a1[1,3]), fmt = "%#.2f"), ")"),
            paste0(sprintf(exp(a1[2,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(a1[2,2] - 1.96*a1[2,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(a1[2,2] + 1.96*a1[2,3]), fmt = "%#.2f"), ")"),
            paste0(sprintf(exp(a1[a1$term %in% c(paste0("t1:", rowname), paste0("t1:", rowname, "TRUE"), paste0("t1TRUE:", rowname), paste0("t1TRUE:", rowname, "TRUE")),2]), fmt = "%#.2f"), " (",
                   sprintf(exp(a1[a1$term %in% c(paste0("t1:", rowname), paste0("t1:", rowname, "TRUE"), paste0("t1TRUE:", rowname), paste0("t1TRUE:", rowname, "TRUE")),2] - 1.96*a1[a1$term %in% c(paste0("t1:", rowname), paste0("t1:", rowname, "TRUE"), paste0("t1TRUE:", rowname), paste0("t1TRUE:", rowname, "TRUE")),3]), fmt = "%#.2f"), "-",
                   sprintf(exp(a1[a1$term %in% c(paste0("t1:", rowname), paste0("t1:", rowname, "TRUE"), paste0("t1TRUE:", rowname), paste0("t1TRUE:", rowname, "TRUE")),2] + 1.96*a1[a1$term %in% c(paste0("t1:", rowname), paste0("t1:", rowname, "TRUE"), paste0("t1TRUE:", rowname), paste0("t1TRUE:", rowname, "TRUE")),3]), fmt = "%#.2f"), ")")),
    p1 = signif(c(a1[1,5], a1[2,5], a1[a1$term %in% c(paste0("t1:", rowname), paste0("t1:", rowname, "TRUE"), paste0("t1TRUE:", rowname), paste0("t1TRUE:", rowname, "TRUE")),5]), 2),
    
    HR2 = c(paste0(sprintf(exp(b1[1,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(b1[1,2] - 1.96*b1[1,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(b1[1,2] + 1.96*b1[1,3]), fmt = "%#.2f"), ")"),
            paste0(sprintf(exp(b1[2,2]), fmt = "%#.2f"), " (",
                   sprintf(exp(b1[2,2] - 1.96*b1[2,3]), fmt = "%#.2f"), "-",
                   sprintf(exp(b1[2,2] + 1.96*b1[2,3]), fmt = "%#.2f"), ")"),
            paste0(sprintf(exp(b1[b1$term %in% c(paste0("t2:", rowname), paste0("t2:", rowname, "TRUE"), paste0("t2TRUE:", rowname), paste0("t2TRUE:", rowname, "TRUE")),2]), fmt = "%#.2f"), " (",
                   sprintf(exp(b1[b1$term %in% c(paste0("t2:", rowname), paste0("t2:", rowname, "TRUE"), paste0("t2TRUE:", rowname), paste0("t2TRUE:", rowname, "TRUE")),2] - 1.96*b1[b1$term %in% c(paste0("t2:", rowname), paste0("t2:", rowname, "TRUE"), paste0("t2TRUE:", rowname), paste0("t2TRUE:", rowname, "TRUE")),3]), fmt = "%#.2f"), "-",
                   sprintf(exp(b1[b1$term %in% c(paste0("t2:", rowname), paste0("t2:", rowname, "TRUE"), paste0("t2TRUE:", rowname), paste0("t2TRUE:", rowname, "TRUE")),2] + 1.96*b1[b1$term %in% c(paste0("t2:", rowname), paste0("t2:", rowname, "TRUE"), paste0("t2TRUE:", rowname), paste0("t2TRUE:", rowname, "TRUE")),3]), fmt = "%#.2f"), ")")),
    p2 = signif(c(b1[1,5], b1[2,5], b1[b1$term %in% c(paste0("t2:", rowname), paste0("t2:", rowname, "TRUE"), paste0("t2TRUE:", rowname), paste0("t2TRUE:", rowname, "TRUE")),5]), 2),
    
    n = rep(dim(t)[1], 3))
  
  print(Sys.time() - x)
  
  return(table)
  
}

# unadjusted
table.int.unadj = NULL
file = "ST10 unadjusted.csv"
for (pred in c("sex", "dm.prev", "htn.prev" ,"smoking", "ALT.high", "obese", "FIB4.cat")) {
  print(pred)
  table.int.unadj = rbind(table.int.unadj, 
                          hr.table.int(t = data,
                                       pred = pred,
                                       rowname = pred))
  write_csv(table.int.unadj, file)
}

# adjusted
table.int.adj = NULL
file = "ST10 adjusted.csv"
for (pred in c("sex", "dm.prev", "htn.prev" ,"smoking", "ALT.high", "obese", "FIB4.cat")) {
  print(pred)
  table.int.adj = rbind(table.int.adj, 
                        hr.table.int(t = data,
                                     pred = pred,
                                     cov = "+ age + sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking + statin1",
                                     rowname = pred))
  write_csv(table.int.adj, file)
}





#### Supp. Table 11: Adjusted associations between TM6SF2 genotype and liver magnetic resonance imaging properties

cov = "age+sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10 + dm.prev + htn.prev + alc.high + smoking+statin1+"
name = "adj-statin"
uv.analysis2(data1 = data0,
             cov = cov, name = name, date = date)
a = uv.analysis2(data1 = data0,
                 cov = cov, name = "adj-statin-all", date = date)
a$cat = "adj-statin-all"

b = uv.analysis2(data1 = data0 %>% filter(FIB4 < 1.3),
                 cov = cov, name = "adj-statin-FIB4.low", date = date)
b$cat = "adj-statin-FIB4.low"

c = uv.analysis2(data1 = data0 %>% filter(FIB4 >= 1.3, FIB4 <= 2.67),
                 cov = cov, name = "adj-statin-FIB4.int", date = date)
c$cat = "adj-statin-FIB4.int"

d = uv.analysis2(data1 = data0 %>% filter(FIB4 > 2.67),
                 cov = cov, name = "adj-statin-FIB4.high", date = date)
d$cat = "adj-statin-FIB4.high"

e = uv.analysis2(data1 = data0 %>% filter(DM == 1),
                 cov = cov, name = "adj-statin-DM", date = date)
e$cat = "adj-statin-DM"
table = rbind(a, b, c, d, e)
write_csv(table[, c(8, 1:7)], "ST11.csv")


#### Figure 2: 10 year cumulative incidence of liver-related events or major adverse cardiovascular events ####
cuminc10y.fig = function(data, rowname = "", 
                         pred1 = "tm6sf2", outcome1 = "time.mace", status1 = "status.mace", 
                         pred2 = "tm6sf2", outcome2 = "time.lre", status2 = "status.lre") {
  
  p = cmprsk::cuminc(ftime = data[, outcome1],
                     fstatus = data[, status1],
                     group = data[, pred1])
  
  p0 = as.data.frame(p$"0 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p0 = p0[p0$time1 == max(p0$time1),]
  p0 = p0[dim(p0)[1],]
  p0 = mutate(p0, 
              inc = 100*est,
              ll = 100*est^(exp(-1.96*sqrt(var)/(est*log(est)))),
              ul = 100*est^(exp(1.96*sqrt(var)/(est*log(est)))))
  
  p1 = as.data.frame(p$"1 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p1 = p1[p1$time1 == max(p1$time1),]
  p1 = p1[dim(p1)[1],]
  p1 = mutate(p1, 
              inc = 100*est,
              ll = 100*est^(exp(-1.96*sqrt(var)/(est*log(est)))),
              ul = 100*est^(exp(1.96*sqrt(var)/(est*log(est)))))
  
  p2 = as.data.frame(p$"2 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p2 = p2[p2$time1 == max(p2$time1),]
  p2 = p2[dim(p2)[1],]
  p2 = mutate(p2,
              inc = 100*est,
              ll = 100*est^(exp(-1.96*sqrt(var)/(est*log(est)))),
              ul = 100*est^(exp(1.96*sqrt(var)/(est*log(est)))))
  
  pcomb = NULL
  pcomb$var = rowname
  pcomb = as.data.frame(pcomb)
  
  pcomb$inc.1.0 = p0$inc
  pcomb$ll.1.0 = p0$ll
  pcomb$ul.1.0 = p0$ul
  
  pcomb$inc.1.1 = p1$inc
  pcomb$ll.1.1 = p1$ll
  pcomb$ul.1.1 = p1$ul
  
  pcomb$inc.1.2 = p2$inc
  pcomb$ll.1.2 = p2$ll
  pcomb$ul.1.2 = p2$ul
  
  # p values based on CC vs TT comparison
  data$pred1.mod1 = ifelse(data[,pred1] == 0, 0, ifelse(data[,pred1] == 2, 1, NA))
  p = cmprsk::cuminc(ftime = data[, outcome1],
                     fstatus = data[, status1],
                     group = data[, "pred1.mod1"])
  pcomb$p.1 = ifelse(p$Tests[1,2] < 0.0001, "<0.0001", signif(p$Tests[1,2], 2))
  
  #
  
  p = cmprsk::cuminc(ftime = data[, outcome2],
                     fstatus = data[, status2],
                     group = data[, pred2])
  
  p0 = as.data.frame(p$"0 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p0 = p0[p0$time1 == max(p0$time1),]
  p0 = p0[dim(p0)[1],]
  p0 = mutate(p0, 
              inc = 100*est,
              ll = 100*est^(exp(-1.96*sqrt(var)/(est*log(est)))),
              ul = 100*est^(exp(1.96*sqrt(var)/(est*log(est)))))
  
  p1 = as.data.frame(p$"1 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p1 = p1[p1$time1 == max(p1$time1),]
  p1 = p1[dim(p1)[1],]
  p1 = mutate(p1, 
              inc = 100*est,
              ll = 100*est^(exp(-1.96*sqrt(var)/(est*log(est)))),
              ul = 100*est^(exp(1.96*sqrt(var)/(est*log(est)))))
  
  p2 = as.data.frame(p$"2 1") %>%
    mutate(time1 = time - 365.25*10) %>%
    dplyr::filter(time1 < 0)
  p2 = p2[p2$time1 == max(p2$time1),]
  p2 = p2[dim(p2)[1],]
  p2 = mutate(p2,
              inc = 100*est,
              ll = 100*est^(exp(-1.96*sqrt(var)/(est*log(est)))),
              ul = 100*est^(exp(1.96*sqrt(var)/(est*log(est)))))
  
  #
  
  pcomb$inc.2.0 = p0$inc
  pcomb$ll.2.0 = p0$ll
  pcomb$ul.2.0 = p0$ul
  
  pcomb$inc.2.1 = p1$inc
  pcomb$ll.2.1 = p1$ll
  pcomb$ul.2.1 = p1$ul
  
  pcomb$inc.2.2 = p2$inc
  pcomb$ll.2.2 = p2$ll
  pcomb$ul.2.2 = p2$ul
  
  # p values based on CC vs TT comparison
  data$pred2.mod1 = ifelse(data[,pred2] == 0, 0, ifelse(data[,pred2] == 2, 1, NA))
  p = cmprsk::cuminc(ftime = data[, outcome2],
                     fstatus = data[, status2],
                     group = data[, "pred2.mod1"])
  pcomb$p.2 = ifelse(p$Tests[1,2] < 0.0001, "<0.0001", signif(p$Tests[1,2], 2))
  
  
  #
  
  return(pcomb)
  
}


tab = cuminc10y.fig(data = data, rowname = "All")
tab = rbind(tab, cuminc10y.fig(data = data %>% filter(FIB4 < 1.3), rowname = "FIB4 < 1.3"))
tab = rbind(tab, cuminc10y.fig(data = data %>% filter(FIB4 >= 1.3, FIB4 <= 2.67), rowname = "FIB4 1.3-2.67"))
tab = rbind(tab, cuminc10y.fig(data = data %>% filter(FIB4 > 2.67), rowname = "FIB4 > 2.67"))

d = NULL
d$Outcome = c(rep("MACE", 12), rep("LRE", 12))
d$Group = rep(c(tab$var[1],
                tab$var[2],
                tab$var[3],
                tab$var[4]), 6)
d$Group = factor(d$Group, levels = tab$var)
d$TM6SF2 = rep(c(rep("CC", 4), rep("CT", 4), rep("TT", 4)), 2)
d=as.data.frame(d)

d$cuminc = c(tab$inc.1.0, tab$inc.1.1, tab$inc.1.2,
             tab$inc.2.0, tab$inc.2.1, tab$inc.2.2)
d$cuminc.l = c(tab$ll.1.0, tab$ll.1.1, tab$ll.1.2,
               tab$ll.2.0, tab$ll.2.1, tab$ll.2.2)
d$cuminc.h = c(tab$ul.1.0, tab$ul.1.1, tab$ul.1.2,
               tab$ul.2.0, tab$ul.2.1, tab$ul.2.2)

dodge = position_dodge(0.9)

delta1 = 
  filter(d, Group == "All", Outcome == "LRE", TM6SF2 == "TT")$cuminc - 
  filter(d, Group == "All", Outcome == "LRE", TM6SF2 == "CC")$cuminc
delta2 = 
  filter(d, Group == "All", Outcome == "MACE", TM6SF2 == "TT")$cuminc - 
  filter(d, Group == "All", Outcome == "MACE", TM6SF2 == "CC")$cuminc
ggplot(data = d %>% filter(Group == "All"),
       aes(y=cuminc, x=Outcome, fill = TM6SF2)) +
  geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(aes(ymin = cuminc.l, ymax = cuminc.h), 
                width = 0.2,
                position = dodge) +
  ylab("10 year cumulative incidence (%)") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 5), limits = c(0, 30)) +
  scale_fill_discrete(expression(italic("TM6SF2"))) +
  theme_bw()+
  ggtitle("All patients") + theme(plot.title = element_text(hjust=.5)) +
  annotate("text", label = paste0("p", tab[1, "p.2"]), x = 1, y = 29, size = 3) +
  annotate("text", label = paste0("p=", tab[1, "p.1"]), x = 2, y = 29, size = 3) +
  annotate("text", label = paste0("Delta = ", round(delta1, 1), "%"), x = 1, y = 28, size = 3) +
  annotate("text", label = paste0("Delta = ", round(delta2, 1), "%"), x = 2, y = 28, size = 3)
ggsave("F2 panel A.pdf", width = 3.5, height = 5)

delta1 = 
  filter(d, Group == "FIB4 < 1.3", Outcome == "LRE", TM6SF2 == "TT")$cuminc - 
  filter(d, Group == "FIB4 < 1.3", Outcome == "LRE", TM6SF2 == "CC")$cuminc
delta2 = 
  filter(d, Group == "FIB4 < 1.3", Outcome == "MACE", TM6SF2 == "TT")$cuminc - 
  filter(d, Group == "FIB4 < 1.3", Outcome == "MACE", TM6SF2 == "CC")$cuminc
ggplot(data = d %>% filter(Group == "FIB4 < 1.3"),
       aes(y=cuminc, x=Outcome, fill = TM6SF2)) +
  geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(aes(ymin = cuminc.l, ymax = cuminc.h), 
                width = 0.2,
                position = dodge) +
  ylab("10 year cumulative incidence (%)") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 5), limits = c(0, 30)) +
  scale_fill_discrete(expression(italic("TM6SF2"))) +
  theme_bw()+
  ggtitle("FIB4 < 1.3") + theme(plot.title = element_text(hjust=.5)) +
  annotate("text", label = paste0("p=", tab[2, "p.2"]), x = 1, y = 29, size = 3) +
  annotate("text", label = paste0("p=", tab[2, "p.1"]), x = 2, y = 29, size = 3) +
  #annotate("text", label = paste0("Delta = ", round(delta1, 1), "%"), x = 1, y = 28, size = 3) +
  annotate("text", label = paste0("Delta = ", round(delta2, 1), "%"), x = 2, y = 28, size = 3)
ggsave("F2 panel B.pdf", width = 3.5, height = 5)

delta1 = 
  filter(d, Group == "FIB4 1.3-2.67", Outcome == "LRE", TM6SF2 == "TT")$cuminc - 
  filter(d, Group == "FIB4 1.3-2.67", Outcome == "LRE", TM6SF2 == "CC")$cuminc
delta2 = 
  filter(d, Group == "FIB4 1.3-2.67", Outcome == "MACE", TM6SF2 == "TT")$cuminc - 
  filter(d, Group == "FIB4 1.3-2.67", Outcome == "MACE", TM6SF2 == "CC")$cuminc
ggplot(data = d %>% filter(Group == "FIB4 1.3-2.67"),
       aes(y=cuminc, x=Outcome, fill = TM6SF2)) +
  geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(aes(ymin = cuminc.l, ymax = cuminc.h), 
                width = 0.2,
                position = dodge) +
  ylab("10 year cumulative incidence (%)") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 5), limits = c(0, 30)) +
  scale_fill_discrete(expression(italic("TM6SF2"))) +
  theme_bw()+
  ggtitle("FIB4 1.3-2.67") + theme(plot.title = element_text(hjust=.5)) +
  annotate("text", label = paste0("p=", tab[3, "p.2"]), x = 1, y = 29, size = 3) +
  annotate("text", label = paste0("p=", tab[3, "p.1"]), x = 2, y = 29, size = 3) +
  annotate("text", label = paste0("Delta = ", round(delta1, 1), "%"), x = 1, y = 28, size = 3) +
  annotate("text", label = paste0("Delta = ", round(delta2, 1), "%"), x = 2, y = 28, size = 3)
ggsave("F2 panel C.pdf", width = 3.5, height = 5)

delta1 = 
  filter(d, Group == "FIB4 > 2.67", Outcome == "LRE", TM6SF2 == "TT")$cuminc - 
  filter(d, Group == "FIB4 > 2.67", Outcome == "LRE", TM6SF2 == "CC")$cuminc
delta2 = 
  filter(d, Group == "FIB4 > 2.67", Outcome == "MACE", TM6SF2 == "TT")$cuminc - 
  filter(d, Group == "FIB4 > 2.67", Outcome == "MACE", TM6SF2 == "CC")$cuminc
ggplot(data = d %>% filter(Group == "FIB4 > 2.67"),
       aes(y=cuminc, x=Outcome, fill = TM6SF2)) +
  geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(aes(ymin = cuminc.l, ymax = cuminc.h), 
                width = 0.2,
                position = dodge) +
  ylab("10 year cumulative incidence (%)") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 5), limits = c(0, 30)) +
  scale_fill_discrete(expression(italic("TM6SF2"))) +
  theme_bw()+
  ggtitle("FIB4 > 2.67") + theme(plot.title = element_text(hjust=.5)) +
  annotate("text", label = paste0("p=", tab[4, "p.2"]), x = 1, y = 29, size = 3) +
  annotate("text", label = paste0("p=", tab[4, "p.1"]), x = 2, y = 29, size = 3) +
  annotate("text", label = paste0("Delta = ", round(delta1, 1), "%"), x = 1, y = 28, size = 3)
#annotate("text", label = paste0("Delta = ", round(delta2, 1), "%"), x = 2, y = 28, size = 3)
ggsave("F2 panel D.pdf", width = 3.5, height = 5)

