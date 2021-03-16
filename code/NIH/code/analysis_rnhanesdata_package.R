# rm(list=ls())
## Please read these comments before running the code contained in this R script!
##
## In order for the code to execute, users need to set the object: dir_supplemental (line 133)
## to be the parent directory of the "code" folder contained in the supplemental material.
## For additional details regarding the analyses performed here, please refer to the manuscript.
##
## The code below is divided into the following "chunks"
##  - Section 0: Install and load all necessary packages including the rnhanesdata package
##               which contains the accelerometry data used in our analysis.
##  - Section 1:
##              1a. Load and merge accelerometry, demographic/comorbidity, and mortality data
##              1b. Create new factor variables which we will use in our analysis.
##                  * We collapse most comorbidity data into either "Yes", "No", or msising from NHANES
##                    questionairre data which allows for "don't know" or "refused" as responses.
##                     We assume individuals who respond "don't know" or "missing" do not have that particular condition.
##                  * We collapse adult education into 3 levels: "less than high school", "high school", and "more than highschool"
##                  * We add a "missing" level to alcohol consumption in order to retain individuals with missing data
##                    for this item in our analysis.
##
##  - Section 2: Calculate commonly used accelerometry summary measures:
##                 * TAC: Total activity counts
##                 * TLAC: Total log(1+activity countss)
##                 * WT: Total wear time
##                 * ST: Total sedentary, sleep, or non-wear time
##                 * MVPA: Total time spent in MVPA
##                 * SATP: Transition probability from sedentary, sleep, or non-wear to active. In the manuscript this is referred to
##                         using the subscript sl/nw
##                 * ASTP: Transition probability from active to sedentary, sleep, or non-wear. In the manuscript this is referred to
##                         using the subscript sl/nw
##               Note that in this section these variables are calculated at the day level. After applying exlcusion criteria
##               in Section 3, we average across days within individuals to get one number per measure for each participant.
##               The exclusion criteria involves excluding days which are deemed to have insufficient wear-time (<10 hours)
##
##  - Section 3: Apply exclusion criteria and create a data frame with one row per subject which will be used as a basis for regression
##               analyses. This data frame is called "data_analysis".
##
##               When estimating complex survey generalized linear models, we create a svydesign() object via the survey package which
##               uses "data_analysis". In order to fit models which use both  the "adjusted" or "unadjusted" survey weights, we create two separate
##               svydesign() objects. This is done in Section 4.d when we perform forward selection.
##
##               Once we've subset the data to obtain "data_analysis", we calculate both adjusted and unadjsuted normalized
##               survey weights. These weights are calculated using the reweight_accel() function
##               (see ?reweight_accel for details). The adjusted weights are calculated using age, gender,
##               and ethnicity strata. The "adjusted" normalizaed weights we use for regression analyses are "wtmec4yr_adj_norm"
##
##               Ecxlusion criteria
##                 * Apply age exclusion (i.e. younger than 50, or 85 and over). Note that individuals age 85 and over
##                   at the time of accelerometer wear have NA (missing) values for the variable RIDAGEEX which records
##                   age in months at the time individuals took part in the exam portion of the study.
##                 * Create a table of pairwise missing data on variables that we intend to include in our prediction model to
##                   see the distribution of missing data. In our analysis individuals are only excluded for:
##                      - Missing body mass index (BMI)
##                      - Missing education
##                      - "Bad" accelerometry data
##                              + fewer than 3 days of accelerometry data with at least 10 hours of estimated wear-time
##                              + device calibration flag recorded by NHANES
##                              + data reliability flag recorded by NHANES
##                      - Missing mortality data
##                      - Individuals recorded as "alive" but had fewer than 5 years of follow-up
##
##  - Section 4: Data analysis
##              4a. Perform (unweighted) functional principal component analysis (FPCA) and survey weighted PCA
##              4b. Use backward selection to identify FPCA features associated with 5-year mortality. Find surrogate measures on
##                  the log transformed activity counts that correlate strongly with the features which are associated with 5-year mortality.
##              4c. Perform scalar on function regression (SoFR) where we include individuals' average (log-transformed) activity profiles
##                  as the functional predictor.
##              4d. Use forward selection to evaluate the predictive value of our set of variables identified in 4b/4c as well as the commonly used
##                  accelerometry features calculated in Section 2 of this code.
##
##                  * Note that functional forms are assumed to be linear and no interactions/effect modifications are considered.
##
##
##
##  By default, the assummption is that this R script will be downloaded with the supplemental material which accompanies
##  the manuscript "Organizing and analyzing the activity in NHANES". This supplemental material is a zipped file with 3 folders:
##    - code
##    - figures
##    - tables
##  Assuming users correctly set the variable "dir_supplemental",
##  running this script will save all figure output to the "figures" folder and
##  .tex files which contain the latex version of tables presented in the manuscript.
##
##  Because the data are relatively large and require a non-trivial amount of working memory (RAM)
##  we "clean up" the workspace as we go, clearing items after the relevant figures/tables/regression results are
##  created/printed to the console. This may need to be taken into consideration if there's a particualr area you'd like to investigate in more detail.
##
##  Finally, although most of the code executes fairly quickly, there are two sections that may take a few minutes to run.
##  The first, calculating survey weighted principal components, will not execute by default. This can be changed by switching
##  make_plot_fpca_vs_svypca to TRUE in Section 0.
##  The other chunk that takes time to run is the forward selection procedure. There's a print statement embedded in the code that will
##  report progress on the procedure, but this may take 15-25 minutes to completely finish.


########################################
##                                    ##
##  Section 0: load required packages ##
##                                    ##
########################################


## Check for packages needed to run analyses/install the rnhanesdata package.
## Note: all these packages are available on CRAN and can therefore be downloaded using the install.packages() function,
##       the rnhanesdata package is not on CRAN due to package size
pckgs <- c("tableone","knitr","kableExtra",   ## packages used for creating Table 1
           "devtools",                        ## package used to download R packages stored on GitHub
           "magrittr","dplyr",                ## packages for merging/transforming data
           "survey",                          ## package used for analyzing complex survey data in R
           "mgcv","refund"                    ## packages used for smoothing/functional regression
           )

sapply(pckgs, function(x) if(!require(x,character.only=TRUE,quietly=TRUE)) {
    install.packages(x)
    require(x, character.only=TRUE)
})
rm(list=c("pckgs"))

## Install the rnhanesdata package and dependencies.
## This may take a few minutes because of the size of the data package.
if(!require("rnhanesdata")){
    install_github("andrew-leroux/rnhanesdata")
    require("rnhanesdata")
}

make_plots  <- TRUE  ## change to FALSE if you don't want to create the figures presented in the manuscript
make_tables <- TRUE  ## change to FALSE if you don't want to create the tables presented in the manuscript

## change "make_plot_fpca_vs_svypca" below to TRUE if you want to plot the first 16 survey weighted principal components versus
## the unweighted functional principal components
make_plot_fpca_vs_svypca <- FALSE


dir_supplemental <- "NIH"  ## directory where supplemental material folder is saved
code_path        <- file.path(dir_supplemental,  "code")     ## file path where helper functions and code to create figures are located
figure_path      <- file.path(dir_supplemental,  "figures")  ## file path where figures will be saved
table_path       <- file.path(dir_supplemental, "tables")   ## file path where tables will be saved

## Source a single helper function: calc_weighted_AUC().
## This  function calculates survey weighted AUC given a set of labels and predictions
source(file.path(code_path,"helper_fns.R"))



#######################################
##                                   ##
##  Section 1a: load and merge data  ##
##                                   ##
#######################################

## load the data
data("PAXINTEN_C");data("PAXINTEN_D")
data("Flags_C");data("Flags_D")
data("Mortality_2011_C");data("Mortality_2011_D")
data("Covariate_C");data("Covariate_D")


## re-code activity counts which are considered "non-wear" to be 0
## this doesn't impact much data, most estimated non-wear times correspond to 0 counts anyway
PAXINTEN_C[,paste0("MIN",1:1440)] <- PAXINTEN_C[,paste0("MIN",1:1440)]*Flags_C[,paste0("MIN",1:1440)]
PAXINTEN_D[,paste0("MIN",1:1440)] <- PAXINTEN_D[,paste0("MIN",1:1440)]*Flags_D[,paste0("MIN",1:1440)]



## Merge covariate, mortality, and accelerometry data
## note that both PAXINTEN_* and Covariate_* have a column
## called "SDDSRVYR" indicating which NHANES wave the data is associated with.
## To avoid duplicating this column in the merged data, we add this variable to the "by"
## argument in left_join()
AllAct_C <- left_join(PAXINTEN_C, Mortality_2011_C, by = "SEQN") %>%
        left_join(Covariate_C, by=c("SEQN", "SDDSRVYR"))
AllAct_D <- left_join(PAXINTEN_D, Mortality_2011_D, by = "SEQN") %>%
        left_join(Covariate_D, by=c("SEQN", "SDDSRVYR"))

AllFlags_C <- left_join(Flags_C, Mortality_2011_C, by = "SEQN") %>%
        left_join(Covariate_C, by=c("SEQN", "SDDSRVYR"))
AllFlags_D <- left_join(Flags_D, Mortality_2011_D, by = "SEQN") %>%
        left_join(Covariate_D, by=c("SEQN", "SDDSRVYR"))

## clean up the workspace for memory purposes
rm(list=c(paste0(c("PAXINTEN_", "Covariate_","Mortality_2011_","Flags_"),rep(LETTERS[3:4],each=4))))

## combine data for the two waves
AllAct   <- rbind.data.frame(AllAct_C,AllAct_D)
AllFlags <- rbind.data.frame(AllFlags_C,AllFlags_D)


## clean up the workspace again
rm(list=c("AllAct_C","AllAct_D","AllFlags_C","AllFlags_D"))




##############################################################################
##                                                                          ##
##  Section 1b: create new variables/relevel factor variables for analyses  ##
##                                                                          ##
##############################################################################

## Code year 5 mortality, NAs for individuals with follow up less than 5 years and alive
AllAct$yr5_mort <- AllFlags$yr5_mort <- as.integer(ifelse(AllAct$permth_exm/12 <= 5 & AllAct$mortstat == 1, 1,
                                                          ifelse(AllAct$permth_exm/12 < 5 & AllAct$mortstat == 0, NA, 0))
                                                    )

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
levels(AllAct$EducationAdult) <- levels(AllFlags$EducationAdult) <- list("Less than high school" = c("Less than 9th grade", "9-11th grade"),
                                                                         "High school" = c("High school grad/GED or equivalent"),
                                                                         "More than high school" = c("Some College or AA degree", "College graduate or above"))

## Re-level alcohol consumption to include a level for "missing"
levels(AllAct$DrinkStatus) <- levels(AllFlags$DrinkStatus) <- c(levels(AllAct$DrinkStatus), "Missing alcohol")
AllAct$DrinkStatus[is.na(AllAct$DrinkStatus)] <- AllFlags$DrinkStatus[is.na(AllAct$DrinkStatus)] <- "Missing alcohol"



## Re-order columns so that activity and wear/non-wear flags are the last 1440 columns of our two
## data matrices. This is a personal preference and is absolutely not necessary.
act_cols <- which(colnames(AllAct) %in% paste0("MIN",1:1440))
oth_cols <- which(!colnames(AllAct) %in% paste0("MIN",1:1440))
AllAct   <- AllAct[,c(oth_cols,act_cols)]
AllFlags <- AllFlags[,c(oth_cols,act_cols)]
rm(list=c("act_cols","oth_cols"))







###########################################################
##                                                       ##
##  Section 2: Calcualte common accelerometery features  ##
##                                                       ##
###########################################################

## Assign just the activity and wear/non-wear flag data to matrices.
## This makes computing the features faster but is technically required.
act_mat  <- as.matrix(AllAct[,paste0("MIN",1:1440)])
flag_mat <- as.matrix(AllFlags[,paste0("MIN",1:1440)])

## replace NAs with 0s
## As described in the manuscript, this only affects 501 minutes for 1 day, for one subject
act_mat[is.na(act_mat)]   <- 0
flag_mat[is.na(flag_mat)] <- 0


AllAct$TAC   <- AllFlags$TAC   <- rowSums(act_mat)
AllAct$TLAC  <- AllFlags$TLAC  <- rowSums(log(1+act_mat))
AllAct$WT    <- AllFlags$WT    <- rowSums(flag_mat)
AllAct$ST    <- AllFlags$ST    <- rowSums(act_mat < 100)
AllAct$MVPA  <- AllFlags$MVPA  <- rowSums(act_mat >= 2020)

## calculate fragmentation measures
bout_mat <- apply(act_mat >= 100, 1, function(x){
    mat <- rle(x)
    sed <- mat$lengths[which(mat$values == FALSE)]
    act <- mat$length[mat$values == TRUE]

    sed <- ifelse(length(sed) == 0, NA, mean(sed))
    act <- ifelse(length(act) == 0, NA, mean(act))
    c(sed,act)
})

AllAct$SBout <- AllFlags$SBout <- bout_mat[1,]
AllAct$ABout <- AllFlags$ABout <- bout_mat[2,]
AllAct$SATP  <- AllFlags$SATP  <- 1/AllAct$SBout
AllAct$ASTP  <- AllFlags$ASTP  <- 1/AllAct$ABout
rm(list=c("act_mat","flag_mat","bout_mat"))



###########################################
##                                       ##
##  Section 3: Apply exclusion criteria  ##
##                                       ##
###########################################


## make dataframe with one row per individual to create table 1.
## Remove columns associated with activity to avoid any confusion.
table_dat <- AllAct[!duplicated(AllAct$SEQN),-which(colnames(AllAct) %in% c(paste0("MIN",1:1440),
                                                                            "TAC","TLAC","WT","ST","MVPA",
                                                                            "SBout","ABout","SATP","ASTP"))]

## subset based on our age inclusion/exclusion criteria
## note that individuals age 85 and over are coded as NA
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
## Education (6) and BMI (35).
criteria_vec <- c("(is.na(table_dat$BMI_cat))",         # missing BMI
                  "(is.na(table_dat$EducationAdult))",  # missing education
                  "(table_dat$SEQN %in% nms_rm)",       # too few "good" days of accelerometery data
                  "((!table_dat$eligstat %in% 1) | is.na(table_dat$mortstat) | is.na(table_dat$permth_exm) | table_dat$ucod_leading %in% \"004\")", # missing mortality data, or accidental death
                  "(table_dat$mortstat == 0 & table_dat$permth_exm/12 < 5)" # less than 5 years of follow up with no mortality
)

## create matrix of pairwise missing data based on our exclusion criterial
tab_miss <- matrix(NA, ncol=length(criteria_vec), nrow=length(criteria_vec))
for(i in seq_along(criteria_vec)){
    for(j in seq_along(criteria_vec)){
        eval(parse(text=paste0("miss_cur <- which(", criteria_vec[i], "&", criteria_vec[j],")")))
        tab_miss[i,j] <- length(miss_cur)
        rm(list=c("miss_cur"))
    }
}
rownames(tab_miss) <- colnames(tab_miss) <- c("BMI","Education","Bad Accel Data","Mortality","Follow-up")
rm(list=c("i","j"))
## view missing data pattern
tab_miss


## add in column indicating exclusion:
##   Exclude = 1 indicates an individual does not meet our inclusion criteria
##   Exclude = 0 indicates an individual does meet our inclusion criteria
eval(parse(text=paste0("table_dat$Exclude <- as.integer(", paste0(criteria_vec,collapse="|"), ")")))


## Create our dataset for analysis with one row per subject
## containing only those subjects who meet our inclusion criteria.
data_analysis  <- subset(table_dat, Exclude == 0)
## get adjusted survey weights using the reweight_accel function
data_analysis  <- reweight_accel(data_analysis)

## Get activity/flag data for only those included participants AND who have 3 good days of data.
## Since we've already removed the "bad" days from Act_Analysis and Act_Flags,
## we need only subset based on subject ID now.
Act_Analysis   <- subset(Act_Analysis, SEQN %in% data_analysis$SEQN)
Flags_Analysis <- subset(Flags_Analysis, SEQN %in% data_analysis$SEQN)

## calculate subject specific averages of the accelerometry features
## using only the "good" days of data
act_var_nms <- c("TAC","TLAC","WT","ST","MVPA","SATP","ASTP")
for(i in act_var_nms){
    data_analysis[[i]] <- vapply(data_analysis$SEQN, function(x) mean(Act_Analysis[[i]][Act_Analysis$SEQN==x],na.rm=TRUE), numeric(1))
}

## verify there's no missingness in the rest of our predictors of interest
vars_interest <- c("Age", "Gender", "Race", "EducationAdult", "SmokeCigs", "DrinkStatus", "BMI_cat",
                   "Diabetes","CHF",  "CHD", "Stroke",
                   "Cancer", "MobilityProblem",
                   "permth_exm")
summary(data_analysis[,c(vars_interest,act_var_nms,"mortstat")])

## clean up the workspace
rm(list=c("AllAct","AllFlags","i","criteria_vec","nms_rm","tab_miss"))
gc()



#####################################################
# Table 1: Compare included vs Excluded individuals #
#####################################################
if(make_tables){
    source(file.path(code_path, "create_table1.R"))
}





###################################
##                               ##
##  Section 4: Data Analysis     ##
##                               ##
###################################


#############################
##                         ##
##  Section 4.a: PCA/fPCA  ##
##                         ##
#############################

## get activity in matrix format, log transform, add 1 because of 0 counts
## There is one individual with 501 minutes recorded as NA. These missing data occur on the last day they wore the device
## for the last 501 minutes of the day. We impute these missing data with 0.
Act <- as.matrix(log(1+Act_Analysis[,paste0("MIN",1:1440)]))
Act[is.na(Act)] <- 0

start <- Sys.time()
fpca_fit <- fpca.face(Act,knots=50)
(end <-  Sys.time()-start)
rm(list=c("start","end"))

## exmine the proportion of variability explained by the first 6 principal components
c(cumsum(fpca_fit$evalues)/sum(fpca_fit$evalues))[1:6]


##calculate the mean and sd of the first 6 PC scores for all participants
for(k in 1:6){
    data_analysis[[paste0("mi",k)]] <- vapply(data_analysis$SEQN, function(x) mean(fpca_fit$scores[Act_Analysis$SEQN==x,k]), numeric(1))
    data_analysis[[paste0("si",k)]] <- vapply(data_analysis$SEQN, function(x) sd(fpca_fit$scores[Act_Analysis$SEQN==x,k]), numeric(1))
}


##########################################
## Figure: first 6 principal components ##
##########################################

if(make_plots){
    source(file.path(code_path,"create_principal_components_plot.R"))
    if(make_plot_fpca_vs_svypca){
        source(file.path(code_path,"create_svyprincipal_components_plot.R"))
    }
}
rm(list=c("make_plot_fpca_vs_svypca"))



############################################################################
##                                                                        ##
##  Section 4.b: Backward selection to identify features of first 6 fPCs  ##
##               associated with 5-year mortality                         ##
##                                                                        ##
############################################################################

## Create a svydesign() object for
## estimating complex survey generalized linear models.
## Here we use the adjusted (re-weighted) 4-year normalized survey weights.
data_analysis_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA,
                               weights = ~wtmec4yr_adj_norm, data = data_analysis, nest = TRUE)

ind_vars  <- c("mi1", "mi2","mi3", "mi4","mi5","mi6",
               "si1","si2","si3","si4","si5", "si6")
inc_vars  <- c("Age", "SmokeCigs", "DrinkStatus", "BMI_cat",
               "Diabetes","CHF",  "CHD", "Stroke",
               "Cancer", "MobilityProblem")
exc_vars <- ind_vars
aic_vec  <- var_vec <- model_vec <- rep(NA, length(ind_vars))
for(i in 1:length(ind_vars)){
    aic_ij  <- rep(NA,length(exc_vars))

    for(k in 1:length(exc_vars)){
        form    <- paste0(c(inc_vars, exc_vars[-k]), collapse="+")
        fit_tmp <- svyglm(as.formula(paste("yr5_mort ~", form)), design=data_analysis_svy,family=quasibinomial())

        aic_ij[k] <- AIC(fit_tmp)[2]

        rm(list=c("fit_tmp","form"))
    }

    k_cur         <- which(aic_ij == min(aic_ij))
    model_vec[i]  <- paste0(c(inc_vars, exc_vars[-k_cur]), collapse="+")
    exc_vars      <- exc_vars[-k_cur]
    aic_vec[i]    <- aic_ij[k_cur]
    rm(list=c("k_cur","aic_ij","k"))
}
## get the final model as the first model where AIC increases after removing a variable
## identified mi1, si1, si5, si6
(backward_model <- model_vec[which(diff(aic_vec) > 0) + 1][1])

fit_logistic_pca <- svyglm(as.formula(paste("yr5_mort ~", backward_model)), design=data_analysis_svy,
                           family=quasibinomial())
summary(fit_logistic_pca)
rm(list=c("aic_vec","i","model_vec","exc_vars","ind_vars","var_vec",
          "inc_vars","data_analysis_svy","backward_model","fit_logistic_pca"))




## Having identified candidate predictors, we now attempt to find surrogate measures for these PC
## features on the original data.
## To qualify as a surrogate, the feature calculated on the raw data should be highly correlated with the PC
## feature it is associated with (~|0.75|+).
## The procedure for finding these surrogate measures involved a fair amount of "guess and check" based on the shape
## of the principal components. The "guessing and checking" is not shown here in the interest of conciseness. We discuss this idea in a bit more detail in the manuscript.
##
## Our search here is by no means exhaustive, and we focused primarily on "average" activity during various
## periods of the day indicated by the shape of each principal component.

#################################
## surrogate for m_{i1}/s_{i1} ##
#################################

## Note that the shape of PC1 is predominantly  a shift in average activity profile
## so we use TLAC as our surrogate measure for m_{i1}, and the day-to-day standard deviation
## of TLAC as our surrogate for s_{i1}
data_analysis$sPC1 <- tapply(Act_Analysis$TLAC, Act_Analysis$SEQN, sd)

## compare association between m_{i1}/s_{i1} with our surrogate measures
cor(data_analysis$TLAC,data_analysis$mi1)
cor(data_analysis$sPC1,data_analysis$si1)



##########################
## surrogate for s_{i5} ##
##########################
t5a <- (11*60+1):(15*60)
t5b <- (16*60+1):(19*60)
s5a <- (rowMeans(Act[,t5a]) - rowMeans(Act[,t5b]))
data_analysis$sPC5 <- tapply(s5a, Act_Analysis$SEQN, sd)
cor(data_analysis$sPC5, data_analysis$si5)
rm(list=c("t5a","t5b","s5a"))



##########################
## surrogate for s_{i6} ##
##########################
t6a <- c((8*60+1):(10*60),
         (15*60+1):(17*60),
         (22*60+1):(24*60))
t6b <-c((5*60+1):(7*60),
        (11*60+1):(13*60),
        (18*60+1):(20*60))
s6a <- rowMeans(Act[,t6a])-rowMeans(Act[,t6b])
data_analysis$sPC6 <- tapply(s6a, Act_Analysis$SEQN, sd)
cor(data_analysis$sPC6,data_analysis$si6)
rm(list=c("t6a","t6b","s6a"))






##########################################################
##                                                      ##
##  Section 4.c:  Scalar on function regression (SoFR)  ##
##                                                      ##
##########################################################

## get individual average activity and day-to-day variability profiles
## using log transformed activity smoothed via fPCA
avg_profiles <- matrix(NA, nrow=nrow(data_analysis),ncol=ncol(Act))
uid          <- unique(data_analysis$SEQN)
for(i in 1:length(uid)){
        avg_profiles[i,] <- colMeans(fpca_fit$Yhat[Act_Analysis$SEQN == uid[i],])
}
## center the functions so that our functional regression is identifiable
## add the centered functions to our dataframe using the I() function which retains matrix structure
## within the dataframe
data_analysis$X_cn <- I(avg_profiles - rep(1, nrow(avg_profiles)) %*% t(colMeans(avg_profiles)))
rm(list=c("uid","i"))


## fit the functional linear model using a cyclic p-spline basis
## Note that although this model accounts for survey weights, it ignores the clustering
## structure.
fit_SoFR <- pfr(yr5_mort ~ Age + BMI_cat + SmokeCigs + DrinkStatus + Race + EducationAdult+
                             CHD + Diabetes + CHF + Stroke +MobilityProblem+ Cancer +
                             lf(X_cn, k=30, bs="cp"),
                   data=data_analysis, family=quasibinomial(), weights=data_analysis$wtmec4yr_adj_norm)

## remove all variable names that got assigned to environment, possibly a bug with pfr?
rm(list=ls()[ls()%in%colnames(data_analysis)])


####################################
## Figure: Functional coefficient ##
####################################
if(make_plots){
    source(file.path(code_path,"create_functional_coefficient_plot.R"))
}

## remove the centered data matrix from our data frame
## this is only done to increase computational speed of the cross-validation below
data_analysis <- data_analysis[,-which(colnames(data_analysis) %in% c("X_cn"))]

rm(list=c("fpca_fit","avg_profiles","fit_SoFR"))








#################################################################################
##                                                                             ##
##  Section 4.d:  Forward selection to identify predictive value of variables  ##
##                                                                             ##
#################################################################################




## In this section we perform forward selection using the adjusted (re-weighted) survey weights,
## unadjusted weights, and no survey weights. We identify variables and output associated with these models
## using the naming convention: _adj, _unadj, _wnwgt
## to denote adjusted, unadjusted, and unweighted results, respectively.



## variables to consider as linear predictors of 5-year mortality
ind_vars <- c("Age", "Gender", "Race", "EducationAdult", "SmokeCigs", "DrinkStatus", "BMI_cat",
              "Diabetes","CHF",  "CHD", "Stroke", "Cancer", "MobilityProblem",
              "sPC1", #"sPC5",
              "sPC6", "ST", "WT", "MVPA","TAC", "TLAC", "SATP", "ASTP")

## examine distrubtion of predictor variables
summary(data_analysis[,ind_vars])

## standardize continuous predictors other than age
## so coefficients are interpretable as
## "a one standard deviation change in ... is associated with ... (in/de)crease in log odds of 5-year mortality"
vars_std <- c("sPC1", "sPC5", "sPC6", "ST", "WT", "MVPA","TAC", "TLAC", "SATP", "ASTP")
for(i in vars_std){
    data_analysis[[i]] <- data_analysis[[i]]/sd(data_analysis[[i]])
}
rm(list=c("vars_std","i"))

## Create svydesign() objects associated with the adjusted and unadjusted
## survey weights. These are the basis for complex survey glm regressions.
data_analysis_svy_adj   <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA,
                                     weights = ~wtmec4yr_adj_norm,
                                     data = data_analysis, nest = TRUE)
data_analysis_svy_unadj <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA,
                                     weights = ~wtmec4yr_unadj_norm,
                                     data = data_analysis, nest = TRUE)


## Create empty vector for incrementing included variables in each of the
## three weighting procedures
inc_vars_adj <- inc_vars_unadj <- inc_vars_unwgt <- c()

## Create empty data frames for storing results for each of the three weighting
## procedures. Each of these data frames contains 4 columns:
##    - Variable: The variable selected in order of inclusion via the forward selection procedure
##    - Cross-Validated AUC: k-fold cross-valided AUC associated with the forward selection procedure after each variable is selected into the model
##    - AIC: Complex survey AIC associated with the forward selection procedure after each variable is selected into the model
##    - AUC: In-sample AUC associated with the forward selection procedure after each variable is selected into the model
auc_mat_adj  <- auc_mat_unadj  <- auc_mat_unwgt  <-
                data.frame("Variable" = rep(NA_character_,length(ind_vars)),
                           "Cross-Validated AUC" = rep(NA_real_,length(ind_vars)),
                           "AIC" = rep(NA_real_,length(ind_vars)),
                           "AUC" = rep(NA_real_,length(ind_vars)),
                           stringsAsFactors = FALSE)


## set the seed so cross-validation results are reproducible
set.seed(1244)
## get the training and testing datasets for 10-fold cross validation
n_folds <- 10
## split the data to have an (approximately) equal number of alive/died in each training/test dataset
inx_id_alive <- which(data_analysis$yr5_mort==0)
inx_id_died  <- which(data_analysis$yr5_mort==1)
nid_alive    <- length(inx_id_alive)
nid_died     <- length(inx_id_died)

inx_ls_alive <- split(sample(inx_id_alive, size=nid_alive, replace=FALSE),
                      rep(1:n_folds,each=ceiling(nid_alive/n_folds))[1:nid_alive])
inx_ls_died <- split(sample(inx_id_died, size=nid_died, replace=FALSE),
                     rep(1:n_folds,each=ceiling(nid_died/n_folds))[1:nid_died])
inx_ls <- lapply(1:n_folds, function(x) c(inx_ls_alive[[x]], inx_ls_died[[x]]))
rm(list=c("inx_id_alive","inx_id_died","nid_alive","nid_died","inx_ls_alive","inx_ls_died"))


## function which will get cut-points to use in AUC calculation
## returns all unique predicted values of the model for a given
get_ctpts <- function(model, type='link',...) c(Inf,sort(unique(predict(model,type=type,...)),decreasing=TRUE))

## Vector to multiply vector of k-fold cross validated weights by to
## get an estimate of cross-valided AUC.
## This accounts for potentially unequal sizes in the test datasets
k_id <- vapply(inx_ls, length, numeric(1))/nrow(data_analysis)

for(i in 1:length(ind_vars)){
        ## get a vector containing all variables not already included in the forward prediction
        ## separately for each of the three weighting procedures.
        exc_vars_adj   <- setdiff(ind_vars, inc_vars_adj)
        exc_vars_unadj <- setdiff(ind_vars, inc_vars_unadj)
        exc_vars_unwgt <- setdiff(ind_vars, inc_vars_unwgt)

        ## create empty matrices/vectors to store AUC/AIC results
        ##   auc_ijk_... Store the cross-valided AUC results
        ##   aic_ij_... Store the complex survey AIC for each of the variables not already included in the model
        ##   auc_ij_... Store the full sample AUC for each of the variables not already included in the model.
        auc_ijk_adj <- auc_ijk_unadj <- auc_ijk_unwgt <- matrix(NA, nrow=length(exc_vars_adj), ncol=n_folds)
        aic_ij_adj  <- aic_ij_unadj <- aic_ij_unwgt  <-
            auc_ij_full_adj <- auc_ij_full_unadj <- auc_ij_full_unwgt <- rep(NA,length(exc_vars_adj))


        ## loop over all variables not currently in the model to calcualte their
        ## improvment to AUC/AIC.
        for(j in 1:length(exc_vars_adj)){
                ## current variable to consider
                var_adj   <- exc_vars_adj[j]
                var_unadj <- exc_vars_unadj[j]
                var_unwgt <- exc_vars_unwgt[j]

                ## create formulas for current regression using
                ## all variables previously identified aand the current variable under consideration
                form_adj   <- paste(c(inc_vars_adj, var_adj), collapse=" +")
                form_unadj <- paste(c(inc_vars_unadj, var_unadj), collapse=" +")
                form_unwgt <- paste(c(inc_vars_unwgt, var_unwgt), collapse=" +")



                ## fit the model to the full data for calculating AIC
                fit_adj   <- svyglm(as.formula(paste("yr5_mort ~", form_adj)),
                                    design=data_analysis_svy_adj, family=quasibinomial())
                fit_unadj <- svyglm(as.formula(paste("yr5_mort ~", form_unadj)),
                                    design=data_analysis_svy_unadj, family=quasibinomial())
                fit_unwgt <- glm(as.formula(paste("yr5_mort ~", form_unwgt)), data=data_analysis, family=binomial())


                ## get AIC
                ## once we have exceeded the design degrees of freedom, we can get
                ## point estimates for our coefficients, but can no longer calculate AIC
                if(fit_adj$df.residual > 0){
                    aic_ij_adj[j]  <- AIC(fit_adj)[2]
                    aic_ij_unadj[j]  <- AIC(fit_unadj)[2]
                }

                aic_ij_unwgt[j]  <- AIC(fit_unwgt)


                ## Calculate weighted and unweighted in-sample AUC using the appropriate weights
                auc_ij_full_adj[j]   <- calc_weighted_AUC(response=predict(fit_adj,type='link'),
                                                          cutpts=get_ctpts(fit_adj,type='link'),
                                                          labels=data_analysis$yr5_mort,
                                                          weights=data_analysis$wtmec4yr_adj_norm)
                auc_ij_full_unadj[j] <- calc_weighted_AUC(response=predict(fit_unadj,type='link'),
                                                          cutpts=get_ctpts(fit_unadj,type='link'),
                                                          labels=data_analysis$yr5_mort,
                                                          weights=data_analysis$wtmec4yr_unadj_norm)
                auc_ij_full_unwgt[j] <- calc_weighted_AUC(response=predict(fit_unwgt,type='link'),
                                                          cutpts=get_ctpts(fit_unwgt,type='link'),
                                                          labels=data_analysis$yr5_mort,
                                                          weights=rep(1, nrow(data_analysis)))

                ## get cross-validated AUC
                for(k in 1:n_folds){
                    ## subset test and training data sets
                    SEQN_train <- data_analysis$SEQN[-inx_ls[[k]]]
                    SEQN_test  <- data_analysis$SEQN[inx_ls[[k]]]
                    data_test  <- subset(data_analysis, SEQN %in% SEQN_test)
                    data_train <- subset(data_analysis, SEQN %in% SEQN_train)


                    ## Fit the appropriate models by subsetting the data.
                    ## By subsetting the existing svydesign objects instead of creating new svydesign objects,
                    ## we retain information on the number of PSU/strata in the original study.
                    fit_adj_cv   <- svyglm(as.formula(paste("yr5_mort ~", form_adj)),
                                           design=subset(data_analysis_svy_adj, SEQN %in% SEQN_train),
                                           family=quasibinomial())
                    fit_unadj_cv <- svyglm(as.formula(paste("yr5_mort ~", form_unadj)),
                                           design=subset(data_analysis_svy_unadj, SEQN %in% SEQN_train),
                                           family=quasibinomial())
                    fit_unwgt_cv <- glm(as.formula(paste("yr5_mort ~", form_unwgt)), data=data_train,
                                        family=quasibinomial())


                    ## Calculate weighted and unweighted AUC using the appropriate weights
                    auc_ijk_adj[j,k]   <- calc_weighted_AUC(response=predict(fit_adj_cv, newdata=data_test, type='link'),
                                                            cutpts=get_ctpts(fit_adj_cv,type='link',newdata=data_test),
                                                            labels=data_test$yr5_mort,
                                                            weights=data_test$wtmec4yr_adj_norm)
                    auc_ijk_unadj[j,k] <- calc_weighted_AUC(response=predict(fit_unadj_cv, newdata=data_test, type='link'),
                                                            cutpts=get_ctpts(fit_unadj_cv,type='link',newdata=data_test),
                                                            labels=data_test$yr5_mort,
                                                            weights=data_test$wtmec4yr_unadj_norm)
                    auc_ijk_unwgt[j,k] <- calc_weighted_AUC(response=predict(fit_unwgt_cv, newdata=data_test, type='link'),
                                                            cutpts=get_ctpts(fit_unwgt_cv,type='link',newdata=data_test),
                                                            labels=data_test$yr5_mort,
                                                            weights=rep(1,nrow(data_test)))


                    rm(list=c("data_train","data_test","SEQN_train","SEQN_test",
                              paste0("fit_",c("adj","unadj","unwgt") ,"_cv")))
                }



                print(j)
                rm(list=c(paste0("fit_",c("adj","unadj","unwgt")),
                          paste0("form_",c("adj","unadj","unwgt")),
                          paste0("var_",c("adj","unadj","unwgt")),
                          "k"))
        }

        ## average across the k-folds using a weighted average of the number of individuals in each test set
        auc_ij_adj   <- auc_ijk_adj %*% k_id
        auc_ij_unadj <- auc_ijk_unadj %*% k_id
        auc_ij_unwgt <- auc_ijk_unwgt %*% k_id

        ## combine results for this iteratin
        auc_j_adj   <- data.frame(exc_vars_adj, auc_ij_adj, aic_ij_adj, auc_ij_full_adj, stringsAsFactors = FALSE)
        auc_j_unadj <- data.frame(exc_vars_unadj, auc_ij_unadj, aic_ij_unadj, auc_ij_full_unadj, stringsAsFactors = FALSE)
        auc_j_unwgt <- data.frame(exc_vars_unwgt, auc_ij_unwgt, aic_ij_unwgt, auc_ij_full_unwgt, stringsAsFactors = FALSE)

        ## identify which variable is associated with the best improvement (or least bad decrease) in CV-AUC
        auc_mat_adj[i,]   <- auc_j_adj[which.max(auc_j_adj[,2]),]
        auc_mat_unadj[i,] <- auc_j_unadj[which.max(auc_j_unadj[,2]),]
        auc_mat_unwgt[i,] <- auc_j_unwgt[which.max(auc_j_unwgt[,2]),]

        ## Create matrices which communicate the predictive power of  all univariate regressions.
        ## These are used to create Table 4 in the manuscript.
        if(i == 1){
            auc_mat_1_adj   <- auc_j_adj[rev(order(auc_j_adj[,2])),]
            auc_mat_1_unadj <- auc_j_unadj[rev(order(auc_j_unadj[,2])),]
            auc_mat_1_unwgt <- auc_j_unwgt[rev(order(auc_j_unwgt[,2])),]
        }

        ## Add in our "best" variable to our list of variabels included in the analysis
        inc_vars_adj   <- c(inc_vars_adj, auc_mat_adj[i,1])
        inc_vars_unadj <- c(inc_vars_unadj, auc_mat_unadj[i,1])
        inc_vars_unwgt <- c(inc_vars_unwgt, auc_mat_unwgt[i,1])


        ## at the end of each iteration, print out the current
        ## adjusted survey weight results
        print(paste(i, " Independent variable finished"))
        print(auc_mat_adj[1:i,]);
        rm(list=c(paste0("auc_j_",c("adj","unadj","unwgt")),
                  paste0("auc_ij_",c("adj","unadj","unwgt")),
                  paste0("auc_ijk_",c("adj","unadj","unwgt")),
                  paste0("auc_ij_full_",c("adj","unadj","unwgt")),
                  paste0("aic_ij_",c("adj","unadj","unwgt")),
                  "j"))
}
rm(list=c("inx_ls",
          paste0("inc_vars_",c("adj","unadj","unwgt")),
          paste0("exc_vars_",c("adj","unadj","unwgt")),
          "i","n_folds","k_id","ind_vars"))



## Final model fits using the adjusted survey weights.
## Fit models using the first 7 predictors (where all models agree) as well as
## the number of predictors identified by AIC and AUC criteria.
inx_svy_adj     <- 7
inx_aic_svy_adj <- min(which(diff(auc_mat_adj[,3]) > 0))
inx_auc_svy_adj <- min(which(diff(auc_mat_adj[,2]) < 0))

fit_final <- svyglm(as.formula(paste("yr5_mort ~", paste0(auc_mat_adj[1:inx_svy_adj,1] ,collapse="+") )),
                        design=data_analysis_svy_adj,
                        family=quasibinomial())
fit_final_aic <- svyglm(as.formula(paste("yr5_mort ~", paste0(auc_mat_adj[1:inx_aic_svy_adj,1] ,collapse="+") )),
                        design=data_analysis_svy_adj,
                        family=quasibinomial())
fit_final_auc <- svyglm(as.formula(paste("yr5_mort ~", paste0(auc_mat_adj[1:inx_auc_svy_adj,1] ,collapse="+") )),
                        design=data_analysis_svy_adj,
                        family=quasibinomial())

## look at the coefficients for the first 7 coefficients
## selected using the adjusted survey weight model
summary(fit_final)


####################################################
# Create Tables:                                   #
#   1) single best predictor ordered by AUC        #
#   2) forward selection procedure                 #
#   3) final mode coefficient estimates and 95% CI #
####################################################
if(make_tables){
    source(file.path(code_path, "create_table_auc.R"))

    ## sourcing "create_table_final_regression.R" will create one warning which can be safely ignored
    source(file.path(code_path, "create_table_final_regression.R"))
}




