rm(list=ls())
library(haven)
library(rnhanesdata)
library(dplyr)
#----------------------------------------------------------------------------
# Process the NHANES Physical Activity Data

# This script includes *some* modifications of code from Andrew Leroux:
  # rnhanesdata: NHANES Accelerometry Data Pipeline. https://github.com/andrew-leroux/rnhanesdata.
# Reference: Leroux et al (2019, Stat in Biosciences). https://link.springer.com/article/10.1007/s12561-018-09229-9
# Reference: Smirnova et al. (2019, J Gerontol A Biol Sci Med Sci). https://www.ncbi.nlm.nih.gov/pubmed/31504213
#----------------------------------------------------------------------------
# Aggregate activity counts into intervals of 'k_min' minutes
k_min = 5 
#----------------------------------------------------------------------------
# Import the data:
#install_github('andrew-leroux/rnhanesdata')
sapply(
  c("rnhanesdata", "devtools","magrittr","dplyr", "forcats", "reshape2"                        
  ), function(x) require(x, character.only=TRUE)
)

# Create a (local) temporary directory 
# where lab measurement (cholesterol, blood pressure) data will be downloaded from the CDC website 
# and then loaded into R. These files need to be downloaded separately as 
# the raw files associated with these lab measurements are not included in the rnhanesdata package
dir.create("tmp")
data_path_dir <- file.path(getwd(), "tmp")

## download the lab measurement data for the cohort 2005-2006
# Total Cholesterol: LBXTC
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/TCHOL_D.XPT", 
              destfile = file.path(data_path_dir,"TCHOL_D.XPT"),method = "curl",cacheOK = TRUE)
# HDL Cholesterol: LBDHDD
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/HDL_D.XPT", 
              destfile = file.path(data_path_dir,"HDL_D.XPT"),method = "curl",cacheOK = TRUE)
# Blood Pressure: BPXSY1 , BPXSY2, BPXSY3 and BPXSY4
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/BPX_D.XPT", 
              destfile = file.path(data_path_dir,"BPX_D.XPT"),method = "curl",cacheOK = TRUE)

## load and merge the lab data
lab_data <- process_covar( varnames=c("LBXTC",  "LBDHDD", "BPXSY1","BPXSY2","BPXSY3", "BPXSY4"), 
                                     localpath=data_path_dir)
## remove temporary directory
unlink(data_path_dir, recursive=TRUE)

## Cardiovascular markers: combine waves (Not necessary: we only use wave D)
CVMarkers <- bind_rows(lab_data$Covariate_D)

## remove participants with missing data for ALL non-ID variables
allMiss   <- apply(CVMarkers, 1, function(x) all(is.na(x[-1])))
CVMarkers <- CVMarkers[!allMiss,]

rm(list=c("lab_data","data_path_dir"))

## load the data
data("PAXINTEN_D") # '05-06 minute level PAXINTEN_D data (1 = Sunday)
data("Flags_D")
data("Covariate_D")

## re-code activity counts which are considered "non-wear" to be 0
## this doesn't impact many data points, most estimated non-wear times correspond to 0 counts
PAXINTEN_D[,paste0("MIN",1:1440)] <-PAXINTEN_D[,paste0("MIN",1:1440)]*Flags_D[,paste0("MIN",1:1440)]

## Merge covariate and accelerometry data
## note that both PAXINTEN_* and Covariate_* have a column
## called "SDDSRVYR" indicating which NHANES wave the data is associated with.
## To avoid duplicating this column in the merged data, we add this variable to the "by"
## argument in left_join()
AllAct <- left_join(PAXINTEN_D, Covariate_D, by=c("SEQN", "SDDSRVYR"))
AllFlags <- left_join(Flags_D, Covariate_D, by=c("SEQN", "SDDSRVYR"))

#merge with cardiovascular markers 
AllAct <- left_join(AllAct, CVMarkers, by = "SEQN")
AllFlags <- left_join(AllFlags, CVMarkers, by = "SEQN")
rm(list=c("CVMarkers"))

## Create Age in years using the age at examination (i.e. when participants wore the device)
AllAct$Age <- AllFlags$Age <- AllAct$RIDAGEEX/12

## Re-level comorbidities to assign refused/don't know as not having the condition
## Note that in practice this does not affect many individuals, but it is an assumption we're making.
levels(AllAct$CHD)    <- levels(AllFlags$CHD)    <- list("No" = c("No","Refused","Don't know"), "Yes" = c("Yes"))
levels(AllAct$CHF)    <- levels(AllFlags$CHF)    <- list("No" = c("No","Refused","Don't know"), "Yes" = c("Yes"))
levels(AllAct$Stroke) <- levels(AllFlags$Stroke) <- list("No" = c("No","Refused","Don't know"), "Yes" = c("Yes"))
levels(AllAct$Cancer) <- levels(AllFlags$Cancer) <- list("No" = c("No","Refused","Don't know"), "Yes" = c("Yes"))
levels(AllAct$Diabetes) <- levels(AllFlags$Diabetes) <- list("No" = c("No","Borderline", "Refused","Don't know"), "Yes" = c("Yes"))

## Re-level education to have 3 levels and categorize don't know/refused to be missing
levels(AllAct$EducationAdult) <- levels(AllFlags$EducationAdult) <- 
  list("Less than high school" = c("Less than 9th grade", "9-11th grade"),
       "High school" = c("High school grad/GED or equivalent"),
       "More than high school" = c("Some College or AA degree", "College graduate or above"))

## Re-level alcohol consumption to include a level for "missing"
levels(AllAct$DrinkStatus) <- levels(AllFlags$DrinkStatus) <- c(levels(AllAct$DrinkStatus), "Missing alcohol")
AllAct$DrinkStatus[is.na(AllAct$DrinkStatus)] <- AllFlags$DrinkStatus[is.na(AllAct$DrinkStatus)] <- "Missing alcohol"

# systolic blood pressure calculation 
AllAct$SYS <- AllFlags$SYS <-  round(apply(AllAct[,c("BPXSY1","BPXSY2","BPXSY3", "BPXSY4")],
                                           1,mean, na.rm= TRUE))

## Re-order columns so that activity and wear/non-wear flags are the last 1440 columns of our two
## data matrices. This is a personal preference and is absolutely not necessary.
act_cols <- which(colnames(AllAct) %in% paste0("MIN",1:1440))
oth_cols <- which(!colnames(AllAct) %in% paste0("MIN",1:1440))
AllAct   <- AllAct[,c(oth_cols,act_cols)]
AllFlags <- AllFlags[,c(oth_cols,act_cols)]
rm(list=c("act_cols","oth_cols"))

# And remove these:
rm(PAXINTEN_D, Covariate_D)

## make dataframe with one row per individual to create table 1.
## Remove columns associated with activity to avoid any confusion.
table_dat <- AllAct[!duplicated(AllAct$SEQN), 
                    -which(colnames(AllAct) %in% 
                             c(paste0("MIN",1:1440),"SBout","ABout","SATP","ASTP"))]

## subset based on our age inclusion/exclusion criteria
## note that individuals age 85 and over are coded as NA
#number of individuals excluded due to subset selection
table_dat <- subset(table_dat, !(Age < 50 | is.na(Age)))

## get the SEQN (id variable) associated with individuals with fewer than 3 days accelerometer wear time
## with at least 10 hours OR had their data quality/device calibration flagged by NHANES
keep_inx       <- exclude_accel(AllAct, AllFlags)
Act_Analysis   <- AllAct[keep_inx,]
Flags_Analysis <- AllFlags[keep_inx,]
nms_rm         <- unique(c(Act_Analysis$SEQN[-which(Act_Analysis$SEQN %in% names(table(Act_Analysis$SEQN))[table(Act_Analysis$SEQN)>=3])],
                           setdiff(AllAct$SEQN, Act_Analysis$SEQN))
)
rm(list=c("keep_inx"))


## Additional inclusion/exclusion criteria.
## Aside from mortality or accelerometer weartime, the only missingness is in
## BMI, Education, SYS, total cholesterol, LBXTC and HDL cholesterol, LBDHDD.
criteria_vec <- c("(is.na(table_dat$BMI_cat))",         # missing BMI
                  "(is.na(table_dat$EducationAdult))",  # missing education
                  "(table_dat$SEQN %in% nms_rm)",       # too few "good" days of accelerometery data
                  "(is.na(table_dat$SYS) | (is.na(table_dat$LBXTC)) | (is.na(table_dat$LBDHDD)) )" #missing lab measures
)

## add in column indicating exclusion:
##   Exclude = 1 indicates an individual does not meet our inclusion criteria
##   Exclude = 0 indicates an individual does meet our inclusion criteria
eval(parse(text=paste0("table_dat$Exclude <- as.integer(", paste0(criteria_vec,collapse="|"), ")")))
table_dat$Exclude <- as.numeric((is.na(table_dat$BMI_cat))|(is.na(table_dat$EducationAdult))|(table_dat$SEQN %in% nms_rm)|(is.na(table_dat$SYS) | (is.na(table_dat$LBXTC)) | (is.na(table_dat$LBDHDD)) ))
## Create our dataset for analysis with one row per subject
## containing only those subjects who meet the inclusion criteria.
data_analysis  <- subset(table_dat, Exclude == 0)

# Rename the survey weights needed:
data_analysis = data_analysis %>% rename("svywt" = "WTMEC2YR")

# Remove the extraneous variables from data analysis 
extra_var = c("PAXCAL","PAXSTAT", "WEEKDAY", "SDDSRVYR", #"WTMEC2YR", 
              "SDMVPSU", "SDMVSTRA", "WTINT2YR", "RIDAGEMN", "RIDAGEEX", "RIDAGEYR",
              "Exclude", "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")
data_analysis = data_analysis %>% select(-extra_var)

# Rename some variables:
data_analysis = data_analysis %>% rename(
  "Total Cholesterol" = "LBXTC",
  "HDL Cholesterol" = "LBDHDD",
  "Systolic Blood Pressure" = "SYS",
)

# Combine some categories of these variables and remove those with mobility problems:
data_analysis = data_analysis %>% mutate(
  Race = case_when(
    Race == "White" ~ "White",
    Race == "Black" ~ "Black",
    Race == "Other" ~ "Other",
    (Race == "Mexican American") | (Race == "Other Hispanic") ~ "Hispanic"
  ) %>% fct_relevel("White", "Black", "Hispanic", "Other"),
  
  Education = case_when(
    EducationAdult == "Less than high school" ~ "< HS",
    EducationAdult == "High school" ~ "= HS",
    EducationAdult == "More than high school" ~ "> HS"
  ) %>% fct_relevel("< HS", "= HS", "> HS"),
  
  MobilityProblem = case_when(
    MobilityProblem == "No Difficulty" ~ "No",
    MobilityProblem == "Any Difficulty" ~ "Yes"
  ) %>% fct_relevel("No", "Yes"),
) %>% filter(
  MobilityProblem == "No"
)

## Get activity/flag data for only those included participants AND days.
## Since we've already removed the "bad" days from Act_Analysis and Act_Flags,
## we need only subset based on subject ID now
Act_Analysis   <- subset(Act_Analysis, SEQN %in% data_analysis$SEQN)
Flags_Analysis <- subset(Flags_Analysis, SEQN %in% data_analysis$SEQN)

## Sort the SEQN in the ascending order in both Act_Analysis and Flags_Analysis
Act_Analysis <- Act_Analysis[sort.int(Act_Analysis$SEQN, index.return = T)$ix,]
Flags_Analysis <- Flags_Analysis[sort.int(Flags_Analysis$SEQN, index.return = T)$ix,]

rm(list=c("AllAct","AllFlags","nms_rm"))
#----------------------------------------------------------------------------
# Construct the subject-specific data
#----------------------------------------------------------------------------
# Number of individuals:
n = length(data_analysis$SEQN)

# Minutes in the day:
times = 1:1440

# Aggregate into intervals:
all_times = rep(seq(k_min, length(times), by = k_min), 
                each = k_min)

# Unique observation points:
tau = sort(unique(all_times)); 

# Compute total activity count by time-of-day:
Act = matrix(NA, nrow = n, ncol = length(tau));

# And number of curve observations per weekday/weekend 
n_wkday = n_wkend = numeric(n)

# Loop over individuals:
for(i in 1:n){
  # Indices of the ith individual:
  inds = which(data_analysis$SEQN[i] == Act_Analysis$SEQN)
  
  # Total activity for individual i in each minute (possibly across multiple days):
  Act_i = colSums(Act_Analysis[inds,paste0("MIN",times)], na.rm=TRUE)
  
  # Aggregate activitiy in k_min-minute intervals:
  Act[i,] = sapply(tau, function(tau_j) sum(Act_i[which(all_times == tau_j)]))
  
  # Count number of weekdays and weekends for individual i:
    # (1 = Sun, 2 = Mon, 3 = Tues, 4 = Weds, 5 = Thurs, 6 = Fri, 7 = Sat)
  wkday_i = Act_Analysis[inds,]$WEEKDAY
  n_wkday[i] = sum((wkday_i>=2) & (wkday_i <= 6))
  n_wkend[i] = sum((wkday_i==1) | (wkday_i == 7))
}

# Add weekday/weekend variables to dataset:
data_analysis = data_analysis %>% 
  mutate(n_wkday = n_wkday,
         n_wkend = n_wkend)

# Subset to complete cases:
sub_subj = complete.cases(Act) & 
  complete.cases(data_analysis)

# Predictor data:
seqn = data_analysis[sub_subj,"SEQN"] 
data_analysis = data_analysis[sub_subj,] %>% select(- c("SEQN"))

# Construct the matrix of predictors:
X = model.matrix( ~ 
                    n_wkday + n_wkend + 
                    BMI + 
                    Race + Gender + Age +
                    Education +
                    DrinksPerWeek + 
                    SmokeCigs + 
                    Diabetes + CHF + CHD + Cancer + Stroke + 
                    `HDL Cholesterol` + 
                    `Total Cholesterol` + 
                    `Systolic Blood Pressure`, data = data_analysis)

rownames(X) = seqn
# Actigraphy count data: raw totals
Ytot = Act[sub_subj,]; colnames(Ytot) = paste(tau, 'min')

# Average activity per minute:
n_days = rowSums(X[,c("n_wkday", "n_wkend")]) # number of days 
m = ncol(Ytot) # number of time points
Y = Ytot/(rep(n_days, m)*k_min)

# Survey weights:
svywt = data_analysis[,'svywt']

# Save the X, Y, Ytot, tau, and the survey weights:
save(X, Y, Ytot, tau, svywt,seqn, file = paste('pa_nhanes_', k_min, 'min.RData', sep=''))