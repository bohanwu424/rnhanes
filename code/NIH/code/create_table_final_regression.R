
## use the fact that the AUC criteria selected the most variables
max_vars <- length(coef(fit_final_auc))

## empty matrix to store coefficient estimates
matrix_reg    <- matrix("", ncol=4, nrow=max_vars)

## add in variable names
matrix_reg[,1] <- names(coef(fit_final_auc))

coef_v   <- round(exp(cbind(coef(fit_final), confint(fit_final))),4)
coef_aic <- round(exp(cbind(coef(fit_final_aic), confint(fit_final_aic))),4)
coef_auc <- round(exp(cbind(coef(fit_final_auc), confint(fit_final_auc))),4)


matrix_reg[1:nrow(coef_v),2]   <- apply(coef_v, 1, function(x) sprintf("%5.3f (%5.3f, %5.3f)",x[1],x[2],x[3]))
matrix_reg[1:nrow(coef_aic),3] <- apply(coef_aic, 1, function(x) sprintf("%5.3f (%5.3f, %5.3f)",x[1],x[2],x[3]))
matrix_reg[1:nrow(coef_auc),4] <- apply(coef_auc, 1, function(x) sprintf("%5.3f (%5.3f, %5.3f)",x[1],x[2],x[3]))



multiple_cat_vars <- c("BMI_cat","DrinkStatus","EducationAdult", "SmokeCigs")

## this will generally produce a warning, can safely ignore
inx_multiple_cat <- sapply(multiple_cat_vars, function(x) min(which(grepl(x,matrix_reg[,1]))))
inx_multiple_cat <- sort(inx_multiple_cat[is.finite(inx_multiple_cat)])

for(p in seq_along(inx_multiple_cat)){
    matrix_reg <- rbind(matrix_reg[1:(inx_multiple_cat[p]-1),],
                        c(names(inx_multiple_cat)[p], "","",""),
                        matrix_reg[(inx_multiple_cat[p]):nrow(matrix_reg),])
    if(p < length(inx_multiple_cat)) inx_multiple_cat[p:length(inx_multiple_cat)] <- inx_multiple_cat[p:length(inx_multiple_cat)] + 1
}
rm(list=c("p","inx_multiple_cat","multiple_cat_vars","max_vars","coef_v","coef_aic","coef_auc"))



## replace variable names with meningful names
matrix_reg <- gsub("DrinkStatus([A-Z].*)", "\\1", matrix_reg)
matrix_reg <- gsub("DrinkStatus", "Alcohol consumption", matrix_reg)
matrix_reg <- gsub("BMI_cat([A-Z].*)", "\\1", matrix_reg)
matrix_reg <- gsub("BMI_cat", "Body mass index", matrix_reg)
matrix_reg <- gsub("SmokeCigs([A-Z].*)", "\\1", matrix_reg)
matrix_reg <- gsub("SmokeCigs", "Cigarette smoking", matrix_reg)
matrix_reg <- gsub("EducationAdult([A-Z].*)", "\\1", matrix_reg)
matrix_reg <- gsub("EducationAdult", "Education", matrix_reg)


matrix_reg <- gsub("^ST", "Sedentary, Sleep, or Non-wear", matrix_reg)
matrix_reg <- gsub("WT", "Wear time", matrix_reg)
matrix_reg <- gsub("MobilityProblem[A-Z].*", "Mobility problem", matrix_reg)
matrix_reg <- gsub("SmokeCigs([A-Z].*)", "\\1", matrix_reg)
matrix_reg <- gsub("CHD[A-Z].*", "Coronary heart disease", matrix_reg)
matrix_reg <- gsub("CHF[A-Z].*", "Congestive heart failure", matrix_reg)
matrix_reg <- gsub("Diabetes[A-Z].*", "Diabetes", matrix_reg)
matrix_reg <- gsub("Cancer[A-Z].*", "Cancer", matrix_reg)
matrix_reg <- gsub("Gender([A-Z].*)", "\\1", matrix_reg)
matrix_reg <- gsub("sPC6", "Surrogate for $s_{i6}$", matrix_reg)
matrix_reg <- gsub("sPC5", "Surrogate for $s_{i5}$", matrix_reg)
matrix_reg <- gsub("sPC1", "Surrogate for $s_{i1}$", matrix_reg)
matrix_reg <- gsub("ASTP", "ASTP$_{sl/nw}$", matrix_reg)
matrix_reg <- gsub("SATP", "SATP$_{sl/nw}$", matrix_reg)

colnames(matrix_reg) <- c("Variable", "", "AUC Stopping Criteria", "AIC Stopping Criteria")



cat("-------------------------------------------------------------------------------------------- \n
    Regression Results from Forward Selection (Adjusted survey weights)  \n
    -------------------------------------------------------------------------------------------- \n ")
print(matrix_reg)


write(kable(matrix_reg, format='latex',booktabs=TRUE,escape=FALSE), file=file.path(table_path, "table_final_regression.tex"))

rm(list=c("matrix_reg"))
