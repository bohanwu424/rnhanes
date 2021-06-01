## create factor variable for 5 year mortality
dd <- data_analysis%>% select(c("SEQN","n_wkday","n_wkend" ,"BMI_cat" , "Race" , "Gender" , "Age" ,"Education" ,
                            "DrinksPerWeek" , "SmokeCigs" , "Diabetes" , "CHF" , "CHD" , "Cancer" , "Stroke" , "HDL Cholesterol" , 
                            "Total Cholesterol" , "Systolic Blood Pressure","svywt"))

mort_dd <- Mortality_2011_D[,c("SEQN","permth_exm", "mortstat")] %>% 
          mutate( yr5_mort =  as.integer(ifelse(permth_exm/12 <= 5 & mortstat == 1, 1,
                                        ifelse(permth_exm/12 < 5 & mortstat == 0, NA, 0))))%>%filter(complete.cases(yr5_mort))%>%
  dplyr::select(c("SEQN", "yr5_mort"  ))

td  <- dd %>%
  inner_join(Covariate_D[,c("SEQN","SDMVPSU","SDMVSTRA")],by = "SEQN")%>%
  inner_join(mort_dd[,c("SEQN","yr5_mort")], by = "SEQN") %>%
  mutate(yr5_mort = factor(yr5_mort, levels=c(0,1), labels=c("Alive","Dead")),
         Age = cut(Age,c(50, 60, 70, 80, 85),right=FALSE),
         svywt = svywt/mean(svywt))


#Processing table_data
vars_interest = setdiff(names(td), c("SEQN","SDMVPSU","SDMVSTRA","svywt","yr5_mort" ))

td_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~svywt, data = td, nest = TRUE)


## get the unweighted and weighted Table 1
tabAsMatrix     <- print(CreateTableOne(vars=c(vars_interest), data =td),
                         printToggle=FALSE,noSpaces=TRUE)[,-4]

tabAsMatrix_svy <- print(svyCreateTableOne(vars= c(vars_interest), data =td_svy),
                         printToggle=FALSE,noSpaces=TRUE)[,-4]
tabAsMatrix_both <- cbind(tabAsMatrix,tabAsMatrix_svy)

## re-name some factors to be interpretable
rownames(tabAsMatrix_both) <- gsub("DrinkStatus", "Alcohol consumption", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("SmokeCigs", "Cigarette smoking", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("EducationAdult", "Education", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("BMI_cat", "Body mass index", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("CHD", "Coronary heart disease", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("CHF", "Congestive heart failure", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("mortstat_fac = Dead \\(%\\)", "Mortality at follow-up (% Dead)", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("= Female \\(%\\)", "(% Female)", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("= Yes \\(%\\)", "(% Yes)", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("Chol_HDL", "HDL Cholesterol", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("Chol_tot", "Total Cholesterol", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("Blood_Pressure", "Blood Pressure", rownames(tabAsMatrix_both))
rownames(tabAsMatrix_both) <- gsub("wk", "week", rownames(tabAsMatrix_both))
colnames(tabAsMatrix_both) <- c("Unweighted","Survey-weighted")

## print table 1
cat("----------------------------------------------------------------------------------------------------- \n
    Table 1: Population characteristics of individuals included in the analysis both unweighted and adjusting for survey weights  \n
    ------------------------------------------------------------------------------------------------------ \n ")
kbl(tabAsMatrix_both)
## latex to create Table 1
table1_tex <- kable(tabAsMatrix_both,format="latex",booktabs=TRUE,longtable =TRUE, caption="Population characteristics of individuals included in the analysis") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
write(table1_tex,file=file.path("Tables","table1.tex"))
rm(list=c("tabAsMatrix_svy","tabAsMatrix","td_svy",
          "td","dd","vars_interest","mort_dd",
          "tabAsMatrix_both","table1_tex"))

## see how many individuals were included
