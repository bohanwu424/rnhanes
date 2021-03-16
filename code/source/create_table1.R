
## create factor variable for 5 year mortality
table_dat$mortstat_fac <- factor(table_dat$mortstat, levels=c(0,1), labels=c("Alive","Dead"))
## create age category variable
table_dat$age_cat <- cut(table_dat$Age,c(50, 60, 70, 80, 85),right=FALSE)

#Processing table_data
table_dat = Covariate_D %>% select(c("SEQN", "BMI_cat", "Race", "Gender", 
                       "Diabetes", "CHF", "CHD", "Cancer","Stroke",
                       "EducationAdult" , "MobilityProblem", "DrinkStatus", 
                       "DrinksPerWeek", "SmokeCigs" ))%>%
  inner_join(table_dat %>% 
               select(c("SDMVPSU","SDMVSTRA","WTMEC2YR","SEQN",
                "mortstat_fac","age_cat","Chol_HDL",
                "Chol_tot","Blood_Pressure"  )),by = "SEQN")

vars_interest = setdiff(names(myData), c("SEQN","SDMVPSU","SDMVSTRA","WTMEC2YR","mortstat_fac" ))

## get survey weighted data using the normalized survey weights
table_dat_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~, data = table_dat, nest = TRUE)

## get the unweighted and weighted Table 1
tabAsMatrix     <- print(CreateTableOne(vars=c(vars_interest,"mortstat_fac"), data =table_dat),
                         printToggle=FALSE,noSpaces=TRUE)[,-4]

## get survey weighted data using the normalized survey weights
table_dat <- reweight_accel(table_dat)
tabAsMatrix_svy <- print(svyCreateTableOne(vars=c(vars_interest,"mortstat_fac"), data =table_dat_svy),
                         printToggle=FALSE,noSpaces=TRUE)[,-4]
tabAsMatrix_both <- cbind(tabAsMatrix,tabAsMatrix_svy)

## re-name some factors to be interpretable
rownames(tabAsMatrix_both) <- gsub("DrinkStatus", "Alcohol consumption", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("SmokeCigs", "Cigarette smoking", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("EducationAdult", "Education", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("BMI_cat", "Body mass index", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("CHD", "Coronary heart disease", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("CHF", "Congestive heart failure", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("MobilityProblem = Any Difficulty \\(%\\)", "Mobility problem (% Any Difficulty)", rownames(tabAsMatrix_both))
#rownames(tabAsMatrix_both) <- gsub("permth_exm", "Months mortality follow-up", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("mortstat_fac = Dead \\(%\\)", "Mortality at follow-up (% Dead)", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("= Female \\(%\\)", "(% Female)", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("= Yes \\(%\\)", "(% Yes)", rownames(tabAsMatrix_both))




## print table 1
cat("----------------------------------------------------------------------------------------------------- \n
    Table 1: Comparing Included vs. Excluded individuals both unweighted and adjusting for survey weights  \n
    ------------------------------------------------------------------------------------------------------ \n ")
print(tabAsMatrix_both)

## latex to create Table 1
table1_tex <- kable(tabAsMatrix_both,format="latex",booktabs=TRUE,caption="Population characteristics of individuals excluded from the analysis") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
write(table1_tex,file=file.path("Tables","table1.tex"))
rm(list=c("tabAsMatrix_svy","tabAsMatrix","table_dat_svy","vars_interest","table1_tex"))


kbl(tabAsMatrix_both)




## see how many individuals were excluded because of accelerometer calibration/suspicious accelerometerdata
## Of the 517 individuals removed for "bad" accelerometry data.
## Individuals in the second row, first column of the table below correspond to said individuals.
##   - 239 were removed for data calibration
##   - 278 were removed for too few days of good data
# table(table_dat$Exclude, table_dat$PAXSTAT %in% 1 & table_dat$PAXCAL %in% 1,useNA='always')
