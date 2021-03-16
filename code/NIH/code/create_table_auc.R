## compare AUC for single best predictor using each of the three weighting procedures
n_pred <- nrow(auc_mat_1_adj)

auc_1 <- cbind(1:n_pred,
               auc_mat_1_adj[,1], sprintf("%5.3f", round(auc_mat_1_adj[,2],3)),
               auc_mat_1_unadj[,1], sprintf("%5.3f", round(auc_mat_1_unadj[,2],3)),
               auc_mat_1_unwgt[,1], sprintf("%5.3f", round(auc_mat_1_unwgt[,2],3))
)
## replace variable names with meningful names
auc_1 <- gsub("DrinkStatus", "Alcohol consumption", auc_1)
auc_1 <- gsub("^ST", "Sedentary, Sleep, or Non-wear", auc_1)
auc_1 <- gsub("WT", "Wear time", auc_1)
auc_1 <- gsub("EducationAdult", "Education", auc_1)
auc_1 <- gsub("BMI_cat", "Body mass index", auc_1)
auc_1 <- gsub("MobilityProblem", "Mobility problem", auc_1)
auc_1 <- gsub("SmokeCigs", "Cigarette smoking", auc_1)
auc_1 <- gsub("CHD", "Coronary heart disease", auc_1)
auc_1 <- gsub("CHF", "Congestive heart failure", auc_1)
auc_1 <- gsub("sPC6", "Surrogate for $s_{i6}$", auc_1)
auc_1 <- gsub("sPC5", "Surrogate for $s_{i5}$", auc_1)
auc_1 <- gsub("sPC1", "Surrogate for $s_{i1}$", auc_1)
auc_1 <- gsub("ASTP", "ASTP$_{sl/nw}$", auc_1)
auc_1 <- gsub("SATP", "SATP$_{sl/nw}$", auc_1)

cat("-------------------------------------------------------------------------------------------- \n
    Single predictor Cross-validated AUC \n
    -------------------------------------------------------------------------------------------- \n ")
print(auc_1)


## get number of variables which would have been selected by either AUC or AIC
## for each of the three different weighting mechanisms
nvars_adj_auc <- min(which(diff(auc_mat_adj[,2]) < 0))
nvars_adj_aic <- min(which(diff(auc_mat_adj[,3]) > 0))

nvars_unadj_auc <- min(which(diff(auc_mat_unadj[,2]) < 0))
nvars_unadj_aic <- min(which(diff(auc_mat_unadj[,3]) > 0))

nvars_unwgt_auc <- min(which(diff(auc_mat_unwgt[,2]) < 0))
nvars_unwgt_aic <- min(which(diff(auc_mat_unwgt[,3]) > 0))


## compare AUC and forward selection for each of the three weighting procedures
auc_full <- cbind(auc_mat_adj[,1], sprintf("%6.3f", round(auc_mat_adj[,2],3)),
                  auc_mat_unadj[,1], sprintf("%6.3f", round(auc_mat_unadj[,2],3)),
                  auc_mat_unwgt[,1], sprintf("%6.3f", round(auc_mat_unwgt[,2],3))
)
## replace variable names with meningful names
auc_full <- gsub("DrinkStatus", "Alcohol consumption", auc_full)
auc_full <- gsub("^ST", "Sedentary, Sleep, or Non-wear", auc_full)
auc_full <- gsub("WT", "Wear time", auc_full)
auc_full <- gsub("EducationAdult", "Education", auc_full)
auc_full <- gsub("BMI_cat", "Body mass index", auc_full)
auc_full <- gsub("MobilityProblem", "Mobility problem", auc_full)
auc_full <- gsub("SmokeCigs", "Cigarette smoking", auc_full)
auc_full <- gsub("CHD", "Coronary heart disease", auc_full)
auc_full <- gsub("CHF", "Congestive heart failure", auc_full)
auc_full <- gsub("sPC6", "Surrogate for $s_{i6}$", auc_full)
auc_full <- gsub("sPC5", "Surrogate for $s_{i5}$", auc_full)
auc_full <- gsub("sPC1", "Surrogate for $s_{i1}$", auc_full)
auc_full <- gsub("ASTP", "ASTP$_{sl/nw}$", auc_full)
auc_full <- gsub("SATP", "SATP$_{sl/nw}$", auc_full)


cat("--------------------------------------------------------------------------------------------- \n
    Forward Selection AUC \n
    -------------------------------------------------------------------------------------------- \n ")
print(auc_full)





auc_full[nvars_adj_aic,1] <- paste0("\\cellcolor{gray!15} ", auc_full[nvars_adj_aic,1])
auc_full[nvars_unadj_aic,3] <- paste0("\\cellcolor{gray!15} ", auc_full[nvars_unadj_aic,3])
auc_full[nvars_unwgt_aic,5] <- paste0("\\cellcolor{gray!15} ", auc_full[nvars_unwgt_aic,5])

auc_full[nvars_adj_auc,1] <- paste0("\\framebox[4cm]{", auc_full[nvars_adj_auc,1],"}")
auc_full[nvars_unadj_auc,3] <- paste0("\\framebox[4cm]{", auc_full[nvars_unadj_auc,3],"}")
auc_full[nvars_unwgt_auc,5] <- paste0("\\framebox[4cm]{", auc_full[nvars_unwgt_auc,5],"}")


# auc_full[nvars_adj_aic,1:2] <- paste0("\\cellcolor{gray!15}\\bf ", auc_full[nvars_adj_aic,1:2])
# auc_full[nvars_unadj_aic,3:4] <- paste0("\\cellcolor{gray!15}\\bf ", auc_full[nvars_unadj_aic,3:4])
# auc_full[nvars_unwgt_aic,5:6] <- paste0("\\cellcolor{gray!15}\\bf ", auc_full[nvars_unwgt_aic,5:6])
#
# auc_full[nvars_adj_auc,1:2] <- paste0("\\cellcolor{blue!15}\\bf ", auc_full[nvars_adj_auc,1:2])
# auc_full[nvars_unadj_auc,3:4] <- paste0("\\cellcolor{blue!15}\\bf ", auc_full[nvars_unadj_auc,3:4])
# auc_full[nvars_unwgt_auc,5:6] <- paste0("\\cellcolor{blue!15}\\bf ", auc_full[nvars_unwgt_auc,5:6])


## create latex versions of both tables
write(kable(auc_1, format='latex',booktabs=TRUE,escape=FALSE), file=file.path(table_path, "table_auc_single_pred.tex"))
write(kable(cbind(1:n_pred,auc_full), format='latex',booktabs=TRUE,escape=FALSE), file=file.path(table_path, "table_auc_forward_pred.tex"))

rm(list=c("auc_full","auc_1","n_pred",
          paste0("nvars_adj_", c("aic","auc")),
          paste0("nvars_unadj_", c("aic","auc")),
          paste0("nvars_unwgt_", c("aic","auc"))
          )
   )
