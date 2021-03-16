## create factor variable for 5 year mortality
td  = inner_join(table_dat, Covariate_D)
td$mortstat_fac <- factor(td$yr5_mort, levels=c(0,1), labels=c("Alive","Dead"))
## create age category variable
td$age_cat <- cut(td$Age,c(50, 60, 70, 80, 85),right=FALSE)

#Processing table_data
td = td %>% 
  select(c("SEQN", "BMI_cat", "Race", "Gender", "Diabetes", 
           "CHF", "CHD", "Cancer","Stroke","EducationAdult" ,
           "MobilityProblem", "DrinkStatus", "DrinksPerWeek", "permth_exm",
           "SmokeCigs" ,"SDMVPSU","SDMVSTRA","weight_adj","SEQN",
            "mortstat_fac","age_cat","Chol_HDL","Chol_tot","Blood_Pressure" ))

vars_interest = setdiff(names(td), c("SEQN","SDMVPSU","SDMVSTRA","weight_adj","mortstat_fac" ))

td_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight_adj, data = td, nest = TRUE)

## get the unweighted and weighted Table 1
tabAsMatrix     <- print(CreateTableOne(vars=c(vars_interest,"mortstat_fac"), data =td),
                         printToggle=FALSE,noSpaces=TRUE)[,-4]


tabAsMatrix_svy <- print(svyCreateTableOne(vars=c(vars_interest,"mortstat_fac"), data =td_svy),
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
rownames(tabAsMatrix_both) <- gsub("permth_exm", "Months mortality follow-up", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("mortstat_fac = Dead \\(%\\)", "Mortality at follow-up (% Dead)", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("= Female \\(%\\)", "(% Female)", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("= Yes \\(%\\)", "(% Yes)", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("Chol_HDL", "HDL Cholesterol", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("Chol_tot", "Total Cholesterol", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("Blood_Pressure", "Blood Pressure", rownames(tabAsMatrix_both))


## print table 1
cat("----------------------------------------------------------------------------------------------------- \n
    Table 1: Comparing Included vs. Excluded individuals both unweighted and adjusting for survey weights  \n
    ------------------------------------------------------------------------------------------------------ \n ")
kbl(tabAsMatrix_both)
## latex to create Table 1
table1_tex <- kable(tabAsMatrix_both,format="latex",booktabs=TRUE,caption="Population characteristics of individuals excluded from the analysis") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
write(table1_tex,file=file.path("..","Tables","table1.tex"))
rm(list=c("tabAsMatrix_svy","tabAsMatrix","td_svy","td","vars_interest","table1_tex"))

## see how many individuals were included
