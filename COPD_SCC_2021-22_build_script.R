#-----------------------------------------------------------#
# COPD SCC 2022
#-----------------------------------------------------------#




sink() # Just putting this here so that if I run it over again it doesn't create more and more sinks...

filename <- "Z:/Group_work/Alex/Encrypted/Alex/COPD/SCC 2022/Logs/COPD_clinical_2021-2022_cleaning_log_"
filedate <- Sys.Date()

sink(file = paste(filename, filedate, ".txt", sep = ""),
     append = FALSE,
     split = TRUE)

cat("\n START \n") # This means that every time I run it it restarts the document instead of getting an
# unuseable document at the end

sink()

sink(file = paste(filename, filedate, ".txt", sep = ""),
     append = TRUE,
     split = TRUE)



library(dplyr)
# library(readstata13)
# library(xlsx)
source("H:/My R functions/MySummary.R")
source("H:/My R functions/lintestOR.R")
source("H:/My R functions/tidyoutput.R")
# library(janitor)
# library(officer)
# library(flextable)
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(ggplot2)
library(survsup)
# library(epitools)
library(psych)
library(lme4)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(car)
library(extrafont)
loadfonts()
fonts()
library(forcats)

tablex <- function(x, y, z) { x %>% select(!!y, !!z) %>% table(useNA = "ifany") }

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

medTableforadmiss <- function(x, varname) {   
  # x is the dataset, varname is the variable name, val is the value of interest (e.g. males) 
  varname <- as.character(varname)
  
  eng <- x %>% filter(country == "England") %>% dplyr::select(varname)
  EN <- sum(eng, na.rm = TRUE)
  engIQR <- round(quantile(eng[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE), 0)
  eng <- paste(engIQR[2], " (", engIQR[1], " to ", engIQR[3], ")", sep = "")
  
  
  wal <- x %>% filter(country == "Wales") %>% dplyr::select(varname)
  WN <- sum(wal, na.rm = TRUE)
  walIQR <- round(quantile(wal[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE), 0)
  wal <- paste(walIQR[2], " (", walIQR[1], " to ", walIQR[3], ")", sep = "")
  
  
  scot <- x %>% filter(country == "Scotland") %>% dplyr::select(varname)
  SN <- sum(scot, na.rm = TRUE)
  scotIQR <- round(quantile(scot[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE), 0)
  scot <- paste(scotIQR[2], " (", scotIQR[1], " to ", scotIQR[3], ")", sep = "")
  
  
  all <- x %>% dplyr::select(varname)
  AN <- sum(all, na.rm = TRUE)
  allIQR <- round(quantile(all[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE), 0)
  all <- paste(allIQR[2], " (", allIQR[1], " to ", allIQR[3], ")", sep = "")
  
  ret <- matrix(c(eng, scot, wal, all), nrow = 1, ncol = 4)
  
  colnames(ret) <- c(paste("England (N=", EN, ")", sep = ""),
                     paste("Scotland (N=", SN, ")", sep = ""),
                     paste("Wales (N=", WN, ")", sep = ""),
                     paste("All (N=", AN, ")", sep = ""))
  
  
  return(ret)
}


meanSumRound <- function(x, variable, roundno) {
  variable <- as.character(variable)
  varcol <- filter(psychic, vars == variable) %>% 
    dplyr::select(vars, N, mean, sd)
  varcol[ ,3:4] <- format(round(varcol[ ,3:4], roundno), nsmall = roundno)
  colnames(varcol) <- paste(variable, colnames(varcol), sep = "_")
  return(varcol[ , -1])
  
}

mediSumRound <- function(x, variable, roundno) {
  variable <- as.character(variable)
  varcol <- filter(psychic, vars == variable) %>% 
    dplyr::select(vars, N, median, lo.quart, hi.quart)
  # function updated so that it just gives numbers back rounded according to roundno,
  # without making any exceptions for midway points etc
  varcol[ ,3:5] <- sprintf(paste0("%.", roundno, "f"), 
                           round(varcol[ ,3:5], roundno), nsmall = roundno) # otherwise use 'roundno'
  
  colnames(varcol) <- paste(variable, colnames(varcol), sep = "_")
  return(varcol[ , -1])
}


FreqSum <- function(x, varname) {
  
  varname <- as.character(varname)
  gen <- x %>% dplyr::select(!!varname) %>% drop_na()
  var_N <- data.frame(nrow(gen))
  colnames(var_N) <- paste0(varname, "_N")
  
  if(nrow(gen) == 0) {return(var_N)}
  
  else {
    
    gen0 <- as.data.frame(table(gen[[1]]))
    gen1 <- as.data.frame(round(prop.table(table(gen[[1]]))*100, 1), nsmall = 1) %>% 
      dplyr::rename(perc = Freq)
    gen2 <- inner_join(gen0, gen1, by = "Var1")
    gen2$perc <- sprintf("%.1f", gen2$perc)
    # gen.E2$England <- paste(gen.E2$Freq, " (", gen.E2$perc, ")", sep = "")
    # gen.E2 <- select(gen.E2, Var1, England)
    for (i in 1:nrow(gen2)) {
      gen3 <- gen2
      gen3$Var1 <- as.character(gen3$Var1)
      gen3 <- gen3[i, ]
      colnames(gen3) <- c("Var1", paste0(varname, "_", gsub(" ", "_", gen3[1,1]), "_n"),
                          paste0(varname, "_", gsub(" ", "_", gen3[1,1]), "_perc")) 
      var_N <- cbind(var_N, gen3[ ,2:3])
    }
    return(var_N)
    
  }
}



medTable <- function(x, varname) {   
  # x is the dataset, varname is the variable name, val is the value of interest (e.g. males) 
  varname <- as.character(varname)
  
  eng <- x %>% filter(country == "England") %>% dplyr::select(varname)
  EN <- length(eng[!is.na(eng)])
  engIQR <- quantile(eng[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  eng <- paste(engIQR[2], " (", engIQR[1], " to ", engIQR[3], ")", sep = "")
  
  
  wal <- x %>% filter(country == "Wales") %>% dplyr::select(varname)
  WN <- length(wal[!is.na(wal)])
  walIQR <- quantile(wal[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  wal <- paste(walIQR[2], " (", walIQR[1], " to ", walIQR[3], ")", sep = "")
  
  
  scot <- x %>% filter(country == "Scotland") %>% dplyr::select(varname)
  SN <- length(scot[!is.na(scot)])
  scotIQR <- quantile(scot[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  scot <- paste(scotIQR[2], " (", scotIQR[1], " to ", scotIQR[3], ")", sep = "")
  
  
  all <- x %>% dplyr::select(varname)
  AN <- length(all[!is.na(all)])
  allIQR <- quantile(all[[1]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  all <- paste(allIQR[2], " (", allIQR[1], " to ", allIQR[3], ")", sep = "")
  
  ret <- matrix(c(eng, scot, wal, all), nrow = 1, ncol = 4)
  
  colnames(ret) <- c(paste("England (N=", EN, ")", sep = ""),
                     paste("Scotland (N=", SN, ")", sep = ""),
                     paste("Wales (N=", WN, ")", sep = ""),
                     paste("All (N=", AN, ")", sep = ""))
  
  
  return(ret)
}

# And another one that will work for calculatng frequencies:

# Changing this so it's inline with what Sophie wants

myFreqTable <- function(x, varname) {
  
  
  varname <- as.character(varname)
  print(varname)
  gen.E <- x %>% filter(country == "England") %>% dplyr::select(!!varname) %>% drop_na()
  EN <- nrow(gen.E)
  gen.E0 <- as.data.frame(table(gen.E[[1]]))
  gen.E1 <- as.data.frame(round(prop.table(table(gen.E[[1]]))*100, 1), nsmall = 1) %>% rename(perc = Freq)
  gen.E2 <- inner_join(gen.E0, gen.E1, by = "Var1")
  gen.E2$England <- paste(format(gen.E2$Freq, big.mark=",", trim=TRUE), " (", # N
                          trimws(format(round(gen.E2$perc, 1), nsmall = 1)), "%)", sep = "") # %
  gen.E2 <- select(gen.E2, Var1, England)
  print(gen.E2)
  
  
  gen.W <- x %>% filter(country == "Wales") %>% dplyr::select(!!varname) %>% drop_na()
  WN <- nrow(gen.W)
  gen.W0 <- as.data.frame(table(gen.W[[1]]))
  gen.W1 <- as.data.frame(round(prop.table(table(gen.W[[1]]))*100, 1), nsmall = 1) %>% rename(perc = Freq)
  gen.W2 <- inner_join(gen.W0, gen.W1, by = "Var1")
  gen.W2$Wales <- paste(format(gen.W2$Freq, big.mark=",", trim=TRUE), " (",
                        trimws(format(round(gen.W2$perc, 1), nsmall = 1)),  "%)", sep = "")
  gen.W2 <- select(gen.W2, Var1, Wales)
  print(gen.W2)
  
  gen.S <- x %>% filter(country == "Scotland") %>% dplyr::select(!!varname) %>% drop_na()
  SN <- nrow(gen.S)
  gen.S0 <- as.data.frame(table(gen.S[[1]]))
  gen.S1 <- as.data.frame(round(prop.table(table(gen.S[[1]]))*100, 1), nsmall = 1) %>% rename(perc = Freq)
  gen.S2 <- inner_join(gen.S0, gen.S1, by = "Var1")
  gen.S2$Scotland <- paste(format(gen.S2$Freq, big.mark=",", trim=TRUE)," (",
                           trimws(format(round(gen.S2$perc, 1), nsmall = 1)),  "%)", sep = "")
  gen.S2 <- select(gen.S2, Var1, Scotland)
  print(gen.S2)
  
  gen.A <- x %>% dplyr::select(!!varname) %>% drop_na()
  AN <- nrow(gen.A)
  gen.A0 <- as.data.frame(table(gen.A[[1]]))
  gen.A1 <- as.data.frame(round(prop.table(table(gen.A[[1]]))*100, 1), nsmall = 1) %>% rename(perc = Freq)
  gen.A2 <- inner_join(gen.A0, gen.A1, by = "Var1")
  gen.A2$All <- paste(format(gen.A2$Freq, big.mark=",", trim=TRUE), " (",
                      trimws(format(round(gen.A2$perc, 1), nsmall = 1)),  "%)", sep = "")
  gen.A2 <- select(gen.A2, Var1, All)
  print(gen.A2)
  
  gen.table <- inner_join(gen.E2, gen.S2, by = "Var1") %>% inner_join(gen.W2, by = "Var1") %>%
    inner_join(gen.A2, by = "Var1")
  colnames(gen.table) <- c(varname, paste("England (N=", format(EN, big.mark=",", trim=TRUE), ")", sep = ""),
                           paste("Scotland (N=", format(SN, big.mark=",", trim=TRUE), ")", sep = ""),
                           paste("Wales (N=", format(WN, big.mark=",", trim=TRUE), ")", sep = ""),
                           paste("All (N=", format(AN, big.mark=",", trim=TRUE), ")", sep = ""))
  
  
  
  # row.names(gen.table) <- gen.table$Var1
  
  return(gen.table)
}




histnorm <- function(g) {
  
  h <- hist(g, breaks = 10, density = 10,
            col = "lightgray", xlab = "Accuracy", main = "Overall") 
  xfit <- seq(min(g, na.rm = TRUE), max(g, na.rm = TRUE), length = 40) 
  yfit <- dnorm(xfit, mean = mean(g, na.rm = TRUE), sd = sd(g, na.rm = TRUE)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  
  plot(h, ylim = c(0, max(yfit)))
  lines(xfit, yfit, col = "black", lwd = 2)
}


nlc <- function(x) {cat(paste("\n", x, "\n", sep = ""))}
CP <- function(x) {write.table(x, "clipboard", sep = "\t", row.names = FALSE)}
CPwithrn <- function(x) {write.table(x, "clipboard", sep = "\t", row.names = TRUE)}




dat <- read.csv("Z:/Group_work/Alex/Encrypted/Secure Enclave_2023-03-29/COPD/SCC_2021-22/data/rawData/COPD-2104-2203-v103-Imperial.csv",
                 header = TRUE, stringsAsFactors = TRUE, na.strings = c("NA", ""))

colnames(dat)


# all.equal(dat$X.Age.At.Admission, dat$X.Age.At.Arrival)

summary(dat$X6.1.Oxygen.Prescribed.During.Admission)
summary(dat$X6.1a.Oxygen.Stipulated.Target.Range)    
summary(dat$X6.1a.1.Other.Oxygen.Stipulated.Target.Range)
summary(dat$X6.2.Oxygen.Administered.During.Admission)


dat <- dat %>% select(study_ID = RowID,
                      patient_ID = PatientID,
                      LSOA = ONSPD_AUG_2021_UK.lsoa11,
                      hosp_code = Org,
                      hosp_name = OrgName,
                      trust_code = Tcode.Now,
                      trust_name = Trust.Now,
                      region = Region,
                      country = Country,
                      arrival_date = X1.2a.Arrival.Date,
                      arrival_time = X1.2b.Arrival.Time,
                      admission_date = X4.1a.Admission.Date,
                      admission_time = X4.1b.Admission.Time,
                      # age = X.Age.At.Arrival, # Just keeps things the same as they were by using admission age
                      age = X.Age.At.Admission,
                      gender = X2.3.Gender,
                      smoke_status = X2.5.Smoking.Status,
                      RSR = X5.1.Respiratory.Review,
                      RSR_date = X5.1a.Date.of.respiratory.review,
                      RSR_time = X5.1b.Time.of.respiratory.review,
                      oxygen_prescribed = X6.1.Oxygen.Prescribed.During.Admission,
                      oxygen_target_range = X6.1a.Oxygen.Stipulated.Target.Range,
                      oxygen_target_range_other = X6.1a.1.Other.Oxygen.Stipulated.Target.Range,
                      oxygen_admin = X6.2.Oxygen.Administered.During.Admission,
                      NIV = X7.1.Acute.treatment.with.NIV,
                      NIV_date = X7.1a.Date.NIV.First.Commenced,
                      NIV_date_NR = X7.1a.1.NIV.Date.Not.Recorded,
                      NIV_time = X7.1b.Time.NIV.first.commenced,
                      NIV_time_NR = X7.1b.1.NIV.Time.Not.Recorded,
                      FEV1_pred_value = X8.1.Most.recently.recorded.FEV1...predicted,
                      FEV1_pred_NR = X8.1.1.FEV1...Predicted.Not.Recorded,
                      FEV1_pred_date = X8.1a.Date.of.last.recorded.FEV1...predicted,
                      FEV1_pred_date_NR = X8.1a.1.FEV1...Predicted.Date.Not.Recorded,
                      FEV1FVC_pred_value = X8.2.Most.recently.recorded.FEV1.FVC.ratio,
                      FEV1FVC_pred_NR = X8.2.1.Recent.FEV1.FVC.Ratio.Not.Recorded,
                      FEV1FVC_pred_date = X8.2a.Date.of.last.recorded.FEV1.FVC.ratio,
                      FEV1FVC_pred_date_NR = X8.2a.1.Last.FEV1.FVC.Ratio.Date.Not.Recorded,
                      life_status = X10.1.Discharge.Life.Status,
                      discharge_date = X10.2.Discharge.Date,
  #                     discharge_time = X6.2b.Discharge.Time, - don't have this for some reason
                      discharge_bundle = X10.3.Discharge.Bundle.Completed,
                      discharge_elements_all = X10.4.Discharge.Elements, # completely useless: drop.
                      DB_inhaler = X10.4.Discharge.Elements...Inhaler.technique.checked,
                      DB_maintenance = X10.4.Discharge.Elements...Medication.issued.classes.reviewed,
                      DB_plan = X10.4.Discharge.Elements...Self.management.plan.provided.or.referred.to.community.team.for.plan,
                      DB_pack = X10.4.Discharge.Elements...Emergency.drug.pack.provided.or.referred.to.community.team.for.pack,
                      DB_unsuitable_for_pack = X10.4.Discharge.Elements...Emergency.drug.pack.not.provided.as.assessed.as.unsuit,
                      DB_oxygen_alert = X10.4.Discharge.Elements...Oxygen.alert.card.provided,
                      DB_smoke = X10.4.Discharge.Elements...Smoking.cessation.drugs.referred.for.behavioral.change.intervention,
                      DB_PR = X10.4.Discharge.Elements...Assessed.for.suitability.for.pulmonary.rehabilitation,
                      DB_FU_72hour = X10.4.Discharge.Elements...Follow.up.requests.at.home.within.72.hours.by.person.or.by.phone,
                      DB_MDT = X10.4.Discharge.Elements...Patient.discussed.at.an.MDT.with.a.community.and.or.primary.care.team,
                      DB_BLF_passport = X10.4.Discharge.Elements...BLF.passport.offered.to.the.patient,
                      DB_none = X10.4.Discharge.Elements...None)
                      # X..Dataset.Version   # useless - drop


colnames(dat)

# dat %>% filter(discharge_bundle == "Patient transferred to another hospital") %>% 
#   select(discharge_elements_all) %>% summary()

# Those transferred can't have anything inputted for the discharge bundle

# Immediately drop all draft records before doing anything else because they will just mess up the dataset
# otherwise.

# dat %>% select(hosp_code, overseas) %>% table()

colnames(dat)

nlc("No. records in original dataset:")
nrow(dat)

# nlc("No. records that are drafts or test hospitals:")
# dat %>% filter(completion_status != "100%" | hosp_code == "YYY") %>% nrow()

nlc("No. records that are drafts or test hospitals:")
dat %>% filter(hosp_code == "YYY") %>% nrow()

# dat <- dat %>% filter(completion_status == "100%" & hosp_code != "YYY")
# nlc("No. records after draft records and test hospitals removed:")
# nrow(dat)

dat <- dat %>% filter(hosp_code != "YYY")
nlc("No. records after draft records and test hospitals removed:")
nrow(dat)





# Need to add in the empty gender factors and smoking status factors, and while we're doing it we might as well
# put it in the correct order

summary(dat$smoke_status)
summary(dat$gender)

dat$gender <- factor(dat$gender, levels = c("Male", "Female", "Transgender", "Not recorded/Preferred not to say"))

summary(dat$smoke_status)

dat$smoke_status <- factor(dat$smoke_status,
                           levels = c("Never smoked", "Ex-smoker", "Current smoker", "Ex-smoker and current vaper",
                                      "Never smoked and current vaper", "Not recorded"))
summary(dat$smoke_status)
# Now we do more building and cleaning:

dat <- dat %>% mutate(study_ID = as.character(study_ID),
                      patient_ID = as.character(patient_ID),
                      LSOA = as.character(LSOA))

# Everything seems to be alright at this point.

# Now would probably be a good point to link to the LSOAs.

# Read in the IMD quintiles


IMDeng <- read.csv("Z:/Group_work/PS_AA/General UK data/IMD/clean_IMD2019_England.csv",
                   header = TRUE, stringsAsFactors = FALSE)
IMDwales <- read.csv("Z:/Group_work/PS_AA/General UK data/IMD/clean_IMD2019_Wales.csv",
                     header = TRUE, stringsAsFactors = FALSE)
IMDscot <- read.csv("Z:/Group_work/PS_AA/General UK data/IMD/clean_IMD2016_Scotland.csv",
                    header = TRUE, stringsAsFactors = FALSE)


# Create the quintiles for the English IMD data

IMDeng$IMD_quintile <- NA

IMDeng$IMD_quintile[IMDeng$IMD_decile == 1] <- 1
IMDeng$IMD_quintile[IMDeng$IMD_decile == 2] <- 1
IMDeng$IMD_quintile[IMDeng$IMD_decile == 3] <- 2
IMDeng$IMD_quintile[IMDeng$IMD_decile == 4] <- 2
IMDeng$IMD_quintile[IMDeng$IMD_decile == 5] <- 3
IMDeng$IMD_quintile[IMDeng$IMD_decile == 6] <- 3
IMDeng$IMD_quintile[IMDeng$IMD_decile == 7] <- 4
IMDeng$IMD_quintile[IMDeng$IMD_decile == 8] <- 4
IMDeng$IMD_quintile[IMDeng$IMD_decile == 9] <- 5
IMDeng$IMD_quintile[IMDeng$IMD_decile == 10] <- 5


IMDeng <- IMDeng %>% select(LSOA = LSOA_code_2011, IMD_quintile_Eng = IMD_quintile)
IMDwales <- IMDwales %>% select(LSOA = LSOA_Code, IMD_quintile_Wal = WIMD_2019_Overall_Quintile)
IMDscot <- IMDscot %>% select(LSOA = LSOA_Code, IMD_quintile_Scot = IMD_quintile)



# Join them together:

dat <- left_join(dat, IMDeng, by = "LSOA")
dat <- left_join(dat, IMDwales, by = "LSOA")
dat <- left_join(dat, IMDscot, by = "LSOA")


dat <- dat %>% mutate(IMD_quintile_Eng = factor(IMD_quintile_Eng),
                      IMD_quintile_Wal = factor(IMD_quintile_Wal),
                      IMD_quintile_Scot = factor(IMD_quintile_Scot))


# create a variable to represent any IMD present

dat <- dat %>% mutate(anyIMD = factor(ifelse(is.na(IMD_quintile_Eng) & is.na(IMD_quintile_Wal) &
                                               is.na(IMD_quintile_Scot), "No IMD", "IMD present")))

# (Or:)
dat %>% filter(is.na(IMD_quintile_Eng) & is.na(IMD_quintile_Wal) & is.na(IMD_quintile_Scot)) %>% nrow()

dat$country <- factor(dat$country, levels = c("England", "Wales"), ordered = FALSE)





# Let's convert all our dates to dates and times to times (as they're currently factors)

dat$arrival_date <- as.Date(dat$arrival_date, "%d/%m/%Y")
dat$admission_date <- as.Date(dat$admission_date, "%d/%m/%Y")
dat$RSR_date <- as.Date(dat$RSR_date, "%d/%m/%Y")
dat$NIV_date <- as.Date(dat$NIV_date, "%d/%m/%Y")
dat$FEV1_pred_date <- as.Date(dat$FEV1_pred_date, "%d/%m/%Y")
dat$FEV1FVC_pred_date <- as.Date(dat$FEV1FVC_pred_date, "%d/%m/%Y")
dat$discharge_date <- as.Date(dat$discharge_date, "%d/%m/%Y")

dat$arrival_time <- as.character(dat$arrival_time)
dat$admission_time <- as.character(dat$admission_time)
dat$RSR_time <- as.character(dat$RSR_time)
dat$NIV_time <- as.character(dat$NIV_time)
dat$FEV1_pred_time <- as.character(dat$FEV1_pred_time)
dat$FEV1FVC_pred_time <- as.character(dat$FEV1FVC_pred_time)

summary(dat$arrival_date)
summary(dat$arrival_time)
summary(dat$NIV_date)
summary(dat$NIV_time)
summary(dat$FEV1_pred_date)
summary(dat$FEV1_pred_time)
summary(dat$discharge_date)
summary(dat$discharge_time)



# sort out the NIV date/time not recorded variables

dat$NIV_time_NR[dat$NIV == "No"] <- NA
dat$NIV_date_NR[dat$NIV == "No"] <- NA

dat %>% select(NIV_date_NR, NIV) %>% table(useNA = "ifany")
dat %>% select(NIV_time_NR, NIV) %>% table(useNA = "ifany")



# We can try converting times to 'difftime' when we try and create the table of when people are given
# particular things.

# One person was admitted in the grey zone of changing time zones! 
# Just give them an extra 15 minutes...

dat$admission_time[dat$admission_date == "2022-03-27" & dat$admission_time == "01:50:00"] <- "02:05:00"


dat <- dat %>% mutate(arrival_DT = as.POSIXct(paste(arrival_date, arrival_time), format="%Y-%m-%d %H:%M:%S"),
                      admission_DT = as.POSIXct(paste(admission_date, admission_time), format="%Y-%m-%d %H:%M:%S"),
                      RSR_DT = as.POSIXct(paste(RSR_date, RSR_time), format="%Y-%m-%d %H:%M:%S"),
                      NIV_DT = as.POSIXct(paste(NIV_date, NIV_time), 
                                                     format="%Y-%m-%d %H:%M:%S"))







head(dat$arrival_DT)
head(dat$NIV_DT)
summary(dat$NIV_DT)

# Life status
summary(dat$life_status)


# Time to admission in hours

dat <- dat %>% mutate(arrival_to_admission_hours = difftime(admission_DT, arrival_DT, units = "hours"))
head(dat$arrival_to_admission_hours)
dat$arrival_to_admission_hours <- as.numeric(dat$arrival_to_admission_hours)

summary(dat$arrival_to_admission_hours)


# Time to NIV in hours

dat <- dat %>% mutate(arrival_to_NIV_hours = difftime(NIV_DT, arrival_DT, units = "hours"))
head(dat$arrival_to_NIV_hours)
dat$arrival_to_NIV_hours <- as.numeric(dat$arrival_to_NIV_hours)

summary(dat$arrival_to_NIV_hours)

dat %>% filter(NIV == "Yes") %>% filter(is.na(NIV_DT)) %>% nrow()
dat %>% filter(NIV == "No") %>% filter(!is.na(NIV_DT)) %>% nrow()

dat %>% filter(NIV == "Yes") %>% filter(is.na(NIV_date) | is.na(NIV_time)) %>% nrow()

dat %>% filter(NIV == "Yes") %>% filter(NIV_date_NR == "Not recorded" | NIV_time_NR == "Not recorded") %>% nrow()

# No descrepancies there.

# create a variable for <2 hours, 2-24 hours, +24 hours, and NR

dat$arrival_to_NIV_cat <- NA
dat$arrival_to_NIV_cat[dat$arrival_to_NIV_hours <2] <- "<2 hours"
dat$arrival_to_NIV_cat[dat$arrival_to_NIV_hours >=2 & dat$arrival_to_NIV_hours < 24] <- "2-24 hours"
dat$arrival_to_NIV_cat[dat$arrival_to_NIV_hours >= 24] <- "24+ hours"
dat$arrival_to_NIV_cat[dat$NIV_date_NR == "Not recorded" | dat$NIV_time_NR == "Not recorded"] <- "NIV given but date or time not recorded"
dat$arrival_to_NIV_cat <- factor(dat$arrival_to_NIV_cat, levels = c("<2 hours", "2-24 hours", "24+ hours", 
                                                                    "NIV given but date or time not recorded"))
summary(dat$arrival_to_NIV_cat)
summary(dat$arrival_to_NIV_hours)
summary(dat$arrival_to_NIV_hours)

dat %>% select(arrival_to_NIV_cat, NIV) %>% table(useNA = "ifany")

dat %>% select(NIV_date_NR, NIV) %>% table(useNA = "ifany")
dat %>% select(NIV_time_NR, NIV) %>% table(useNA = "ifany")


# FEV1 sort out

dat$FEV1_pred <- "No"
dat$FEV1_pred[!is.na(dat$FEV1_pred_value)] <- "Yes"
dat$FEV1_pred <- factor(dat$FEV1_pred, levels = c("No", "Yes"))

summary(dat$FEV1_pred)
dat %>% select(FEV1_pred, FEV1_pred_date_NR) %>% table(useNA = "ifany")
dat$FEV1_pred_date_NR[dat$FEV1_pred == "No"] <- NA
dat %>% select(FEV1_pred, FEV1_pred_date_NR) %>% table(useNA = "ifany")
dat %>% filter(!is.na(FEV1_pred_date)) %>% filter(FEV1_pred_date_NR == "Not recorded") %>% nrow()



# FEV1FVC sort out

dat$FEV1FVC_pred <- "No"
dat$FEV1FVC_pred[!is.na(dat$FEV1FVC_pred_value)] <- "Yes"
dat$FEV1FVC_pred <- factor(dat$FEV1FVC_pred, levels = c("No", "Yes"))

summary(dat$FEV1FVC_pred)
dat %>% select(FEV1FVC_pred, FEV1FVC_pred_date_NR) %>% table(useNA = "ifany")
dat$FEV1FVC_pred_date_NR[dat$FEV1FVC_pred == "No"] <- NA
dat %>% select(FEV1FVC_pred, FEV1FVC_pred_date_NR) %>% table(useNA = "ifany")
dat %>% filter(!is.na(FEV1FVC_pred_date)) %>% filter(FEV1FVC_pred_date_NR == "Not recorded") %>% nrow()



# any spirometry available


dat$spirometry <- "No"
dat$spirometry[!is.na(dat$FEV1_pred_value) | !is.na(dat$FEV1FVC_pred_value)] <- "Yes"
dat$spirometry <- factor(dat$spirometry, levels = c("No", "Yes"))
summary(dat$spirometry)

dat$obstruction <- NA
dat$obstruction[dat$FEV1FVC_pred_value < 0.7] <- "Yes (<0.7)"
dat$obstruction[dat$FEV1FVC_pred_value >= 0.7] <- "No (>=0.7)"
dat$obstruction <- factor(dat$obstruction, levels = c("No (>=0.7)", "Yes (<0.7)"))

colnames(dat)



summary(dat$obstruction)


dat %>% filter(FEV1_pred_date > arrival_date) %>% nrow()
dat %>% filter(FEV1FVC_pred_date > arrival_date) %>% nrow()

# length of stay days

dat <- dat %>% mutate(LOS_days = discharge_date - arrival_date)

head(sort(dat$LOS_days))


# Time to RSR in hours

dat <- dat %>% mutate(admission_to_RSR_hours = difftime(RSR_DT, admission_DT, units = "hours"))
head(dat$admission_to_RSR_hours)
dat$admission_to_RSR_hours <- as.numeric(dat$admission_to_RSR_hours)

summary(dat$discharge_bundle)


# Need to remove those who were transferred from this

dat <- dat %>% mutate(LOS_days = ifelse(life_status == "Died as inpatient", 
                                         NA, LOS_days))

summary(dat$oxygen_admin)
summary(dat$oxygen_prescribed)

table(dat$oxygen_admin, dat$oxygen_prescribed)

table(dat$oxygen_target_range, dat$oxygen_prescribed, useNA = "ifany")

# For oxygen administered, we want the denominator to be oxygen prescribed.

dat$oxygen_admin[dat$oxygen_prescribed == "No" ] <- NA

# # # # 


# day of arrival and day of discharge


dat <- dat %>% mutate(arrival_day_of_week = weekdays(arrival_date, abbreviate = FALSE))

dat$arrival_day_of_week <- ordered(dat$arrival_day_of_week, levels=c("Monday", "Tuesday", "Wednesday", "Thursday", 
                                                   "Friday", "Saturday", "Sunday"))

table(dat$arrival_day_of_week)


dat <- dat %>% mutate(discharge_day_of_week = weekdays(discharge_date, abbreviate = FALSE))

str(dat$discharge_day_of_week)

dat$discharge_day_of_week[dat$discharge_day_of_week %in% c("Monday", "Tuesday", "Wednesday", "Thursday", 
                                                           "Friday")] <- "Weekday"

dat$discharge_day_of_week[dat$discharge_day_of_week %in% c("Saturday", "Sunday")] <- "Weekend" 

dat$discharge_day_of_week <- factor(dat$discharge_day_of_week, levels=c("Weekday", "Weekend"))




# Those who transferred need to be removed from the discharge day of week

# Using ifelse here changes them to factor levels so I'm doing it the old school method

# dat <- dat %>% mutate(discharge_day_of_week = ifelse(discharge_bundle == "Patient transferred to another hospital", 
#                                          NA, discharge_day_of_week))



dat$discharge_day_of_week[dat$life_status == "Died as inpatient"] <- NA

summary(dat$life_status)
summary(dat$discharge_bundle)
summary(dat$discharge_day_of_week)

# time of arrival for the time of arival table.
# basically times are so annoying to work with, so what I've done is given each time an associated date of 
# 1970-01-01 (the common origin date) and then subtracted the time from 1970-01-01 00:00:00 to give minutes
# into the day. Seems to have worked well

dat <- dat %>% mutate(arrival_time_mins_from_midnight = as.POSIXct(paste("1970-01-01", arrival_time))) %>%
  mutate(arrival_time_mins_from_midnight = difftime(arrival_time_mins_from_midnight, 
                                                    "1970-01-01 00:00:00", units = "mins"))
head(dat$arrival_time_mins_from_midnight)

dat %>% select(arrival_time, arrival_time_mins_from_midnight) %>% 
  arrange(arrival_time_mins_from_midnight) %>% head(10)
dat %>% select(arrival_time, arrival_time_mins_from_midnight) %>% 
  arrange(desc(arrival_time_mins_from_midnight)) %>% head(10)

dat <- dat %>% mutate(arrival_time_mins_from_midnight = as.integer(arrival_time_mins_from_midnight))



dat$weekday_admission_daynight <- "Weekday day admission"
dat$weekday_admission_daynight[dat$arrival_time_mins_from_midnight < (8*60)] <- "Weekday night admission"
dat$weekday_admission_daynight[dat$arrival_time_mins_from_midnight >= (18*60)] <- "Weekday night admission"

dat$weekday_admission_daynight[dat$arrival_day_of_week %in% c("Saturday", "Sunday")] <- "Weekend day admission"
dat$weekday_admission_daynight[dat$arrival_day_of_week %in% c("Saturday", "Sunday") & dat$arrival_time_mins_from_midnight < (8*60)] <- "Weekend night admission"
dat$weekday_admission_daynight[dat$arrival_day_of_week %in% c("Saturday", "Sunday") & dat$arrival_time_mins_from_midnight >= (18*60)] <- "Weekend day admission"

dat$weekday_admission_daynight <- factor(dat$weekday_admission_daynight, levels = c("Weekday day admission", "Weekday night admission", 
                                                                                    "Weekend day admission", "Weekend night admission"))
summary(dat$weekday_admission_daynight)



dat <- dat %>% mutate(weekday_admission = factor(ifelse( (arrival_day_of_week == "Friday" & arrival_time_mins_from_midnight >= (17*60)) |
                                                           arrival_day_of_week == "Saturday" |
                                                           arrival_day_of_week == "Sunday" |
                                                           (arrival_day_of_week == "Monday" & arrival_time_mins_from_midnight < (8*60)),
                                                         "Weekend admission", "Weekday admission"),
                                                 levels = c("Weekday admission", "Weekend admission")))

dat$daynight_admission <- "Day time"
dat$daynight_admission[dat$arrival_time_mins_from_midnight < (8*60)] <- "Night time"
dat$daynight_admission[dat$arrival_time_mins_from_midnight >= (17*60)] <- "Night time"

dat$daynight_admission <- factor(dat$daynight_admission, levels = c("Day time", "Night time"))


dat$weekday_daynight_admission <- "Weekday day time"
dat$weekday_daynight_admission[dat$weekday_admission == "Weekday admission" & dat$daynight_admission == "Night time"] <- "Weekday night time"
dat$weekday_daynight_admission[dat$weekday_admission == "Weekend admission" & dat$daynight_admission == "Day time"] <- "Weekend day time"
dat$weekday_daynight_admission[dat$weekday_admission == "Weekend admission" & dat$daynight_admission == "Night time"] <- "Weekend night time"

dat$weekday_daynight_admission <- factor(dat$weekday_daynight_admission, levels = c("Weekday day time", "Weekday night time", 
                                                                                    "Weekend day time", "Weekend night time"))
summary(dat$weekday_daynight_admission)

summary(dat$daynight_admission)

table(dat$daynight_admission, dat$weekday_admission)



# Discharge bundle consistency check

dat %>% select(discharge_bundle, DB_none) %>% table(useNA = "ifany") 
table(dat$DB_none)




dat <- dat %>% mutate_at(.vars = vars(starts_with("DB")),
                         .funs = ~ifelse(life_status == "Died as inpatient", NA, .))


dat$DB_smoke[dat$smoke_status != "Current smoker" | is.na(dat$smoke_status)] <- NA
dat %>% select(smoke_status, DB_smoke) %>% table(useNA = "ifany")


summary(dat$discharge_bundle)

dat$discharge_bundle <- factor(dat$discharge_bundle, levels = c("Yes", "No", "Self discharge"))

dat$discharge_day_of_week[dat$life_status != "Alive"] <- NA
summary(dat$discharge_day_of_week)



# Extra variables:

summary(dat$anyIMD)

# IMD from any country:

dat$IMD_quintile_all <- dat$IMD_quintile_Eng
dat$IMD_quintile_all[is.na(dat$IMD_quintile_all)] <- dat$IMD_quintile_Wal[is.na(dat$IMD_quintile_all)]
dat$IMD_quintile_all[is.na(dat$IMD_quintile_all)] <- dat$IMD_quintile_Scot[is.na(dat$IMD_quintile_all)]



summary(dat$discharge_day_of_week)


# Filter out invalid dates (N = 0)

dat %>% filter(discharge_date < "2021-04-01") %>% nrow()
dat %>% filter(discharge_date > "2022-03-31") %>% nrow()

# Was anyone discharged before they arrived, or received medication before they arrived, or received medication
# after they left? No.


nlc("This many people have an admission date before they arrive:")
dat %>% filter(admission_date - arrival_date < 0) %>% nrow() %>% nlc() # 0

# dat <- dat %>% filter(!(admission_date - arrival_date < 0))

nlc("This many people have an NIV date before they arrived:")
dat %>% filter(NIV_DT - arrival_DT < 0) %>% nrow() %>% nlc() # 0

nlc("This many people have an RSR date before they arrived:")
dat %>% filter(RSR_date - arrival_date < 0) %>% nrow() %>% nlc() # 0

nlc("This many people have a discharge date before they arrived:")
dat %>% filter(discharge_date - arrival_date < 0) %>% nrow() %>% nlc() 

dat <- dat %>% filter(!(discharge_date - arrival_date < 0))

nlc("Removed to leave this many:")
nrow(dat) %>% nlc()



# remove duplicates

nlc("assess whether there are duplicate records. Done based on:
country, trust_code, patient_ID, hosp_code, LSOA, age, gender, arrival_date, discharge_date")

nlc("We have this many duplicated records:")

dat %>% select(country, trust_code, patient_ID, hosp_code, LSOA, age, gender, arrival_date, discharge_date) %>% duplicated() %>% 
  sum() %>% nlc()

nrow(dat)




dat <- dat %>% arrange(patient_ID, arrival_DT)



dat <- dat[!duplicated(select(dat, country, trust_code, patient_ID, hosp_code, LSOA, age, gender, arrival_date, discharge_date)), ]


nlc("These are removed to leave this many records:")
nrow(dat) %>% nlc()


# # and do the arrival times
# 
# dat$arrival8hourtimes <- cut(dat$arrival_time_mins_from_midnight, breaks = seq(-0.5, 1439.5, 480),
#                              labels = paste0("lessthan", seq(8, 24, 8)))
# 
# dat$weekday_8hours <- interaction(dat$arrival_day_of_week, dat$arrival8hourtimes, sep = "_")
# 
# # and the weekday arrivals
# 
# summary(dat$weekday_admission)

dat %>% select(starts_with("DB")) %>% summary()

dat$RSR_24hour <- NA
dat$RSR_24hour[dat$admission_to_RSR_hours < 24] <- "<24 hours"
dat$RSR_24hour[dat$admission_to_RSR_hours >= 24] <- ">=24 hours"
dat$RSR_24hour <- factor(dat$RSR_24hour)
summary(dat$RSR_24hour)



# create a BPT variable

summary(dat$discharge_bundle)

dat %>% select(starts_with("DB_")) %>% colnames()

# create the combo DB variables

dat$DB_plan_or_pack_or_unsuitable <- 0
dat$DB_plan_or_pack_or_unsuitable[dat$DB_plan == 1 | dat$DB_pack == 1 | dat$DB_unsuitable_for_pack == 1] <- 1
dat$DB_plan_or_pack_or_unsuitable[is.na(dat$DB_plan)] <- NA


dat <- dat %>% mutate(BPT_all = factor(ifelse(DB_inhaler == 1 & DB_maintenance == 1 & 
                                                (DB_plan == 1 | DB_pack == 1 | DB_unsuitable_for_pack == 1) &
                                                (DB_smoke == 1 | is.na(DB_smoke)) & DB_PR == 1 & 
                                                DB_FU_72hour & discharge_bundle == "Yes" & 
                                                RSR == "Yes", "Achieved", "Not achieved")))

dat$BPT_all[dat$life_status == "Died"] <- NA

table(dat$BPT_all) %>% prop.table()

dat %>% select(DB_smoke, smoke_status) %>% table(useNA = "ifany")
                      
table(dat$DB_FU_72hour)
table(dat$DB_PR)
table(dat$DB_plan)





saveRDS(dat, "Z:/Group_work/Alex/Encrypted/Alex/COPD/SCC 2022/data/tidyData/COPD_SCC_2021-2022_clean_data_2022-09-22.RDS")


