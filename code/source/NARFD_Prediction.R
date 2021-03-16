#---------------------------------------
#NARFD Prediction
load(here::here("data", "NARFD_results.RData"))

## Make Score Dataframe
scores_m = reshape2::dcast(narfd_results$scores, .id ~ prototype,value.var = "Value") %>% column_to_rownames(".id")%>% apply(MARGIN = 2, FUN = scale)
colnames(scores_m) = paste0("score",seq(6))

##Rename Some Variables
X2 = as.data.frame(X)
names(X2)[2:7] = paste0("wk_day_",seq(2,7))
names(X2)[-seq(1,ncol(X2) - 3)] = c("Chol_HDL", "Chol_tot","Blood_Pressure")
names(X2)[14:15] = c("HS", "Above_HS")
names(X2)[1] = "Intercept"

#Include SEQN
X2[,"SEQN"] = seqn

##Merging X with Scores,SEQN,Covariates,Mortaility, Weighted Adj
table_dat <- X_df <- as.data.frame(cbind(X2,scores_m)) %>%
  inner_join(Covariate_D[,c("SEQN","SDMVPSU","SDMVSTRA","WTMEC2YR")],by = "SEQN")%>% 
  inner_join(Mortality_2011_D, by = "SEQN")%>%
  mutate(weight_adj = WTMEC2YR/mean(WTMEC2YR)) %>%
  mutate( yr5_mort =  as.integer(ifelse(permth_exm/12 <= 5 & mortstat == 1, 1,
                                        ifelse(permth_exm/12 < 5 & mortstat == 0, NA, 0))))%>%
  filter(complete.cases(yr5_mort)) 

## Create a svydesign() object
X_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight_adj, data = X_df, nest = TRUE)

#---------------------------------------------------------------
#AIC Selection
ind_vars  <- paste0("score",seq(6))
inc_vars  <- colnames(X2)
exc_vars <- ind_vars
# aic_vec  <- var_vec <- model_vec <- rep(NA, length(ind_vars))
# for(i in 1:length(ind_vars)){
#   aic_ij  <- rep(NA,length(exc_vars))
# 
#   for(k in 1:length(exc_vars)){
#     form    <- paste0(c(inc_vars, exc_vars[-k]), collapse="+")
#     fit_tmp <- svyglm(as.formula(paste("yr5_mort ~", form)), design=X_svy,family=binomial())
#     aic_ij[k] <- fit_tmp$aic
#     rm(list=c("fit_tmp","form"))
#   }
# 
#   k_cur         <- which(aic_ij == min(aic_ij))
#   model_vec[i]  <- paste0(c(inc_vars, exc_vars[-k_cur]), collapse="+")
#   exc_vars      <- exc_vars[-k_cur]
#   aic_vec[i]    <- aic_ij[k_cur]
#   rm(list=c("k_cur","aic_ij","k"))
# }
# 
# ## get the final model as the first model where AIC increases after removing a variable
# (backward_model <- model_vec[which(diff(aic_vec) > 0) + 1][1])
# if(is.na(backward_model)) backward_model = model_vec[1][1]

backward_model <- paste0(c(inc_vars, exc_vars), collapse="+")
rm(list = c("X2","X_df", "score_m","scores_m","samples","ind_vars","inc_vars","exc_vars"))

#---------------------------------------------
## estimating complex survey generalized linear models.
fit_narfd <- svyglm(as.formula(paste("yr5_mort ~", backward_model)), design=X_svy,  family=binomial())
#--------------------------------------------------
##ROC Curves, AUC
pred_ROCR = prediction(fit_narfd$fitted.values, fit_narfd$y)

roc_ROCR = performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
AUC = performance(pred_ROCR, measure = "auc")
AUC = AUC@y.values[[1]]

confusionMatrix(as.factor(1*fit_narfd$fitted.values>0.5), as.factor(fit_narfd$y==1))

rm(list = c("fit_narfd", "pred_ROCR","X_svy","backward_model"))
