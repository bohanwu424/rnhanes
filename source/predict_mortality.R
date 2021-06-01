#---------------------------------------------------------
source(here::here("source", "narfd_create_mat.R")) #create score mat
#---------------------------------------------------------
#Create Train/Test Split
##Rename Some Variables
X2 = as.data.frame(X)
names(X2) = gsub("Education= HS","HS",names(X2) )
names(X2) = gsub("Education> HS","Above_HS",names(X2) )
names(X2) = gsub("`HDL Cholesterol`","Chol_HDL", names(X2))
names(X2) = gsub("`Total Cholesterol`","Chol_tot", names(X2))
names(X2) = gsub("`Systolic Blood Pressure`","Blood_Pressure", names(X2))
names(X2)[1] = "intercept"
#Include SEQN
X2[,"SEQN"] = seqn
#Merging X with Scores,SEQN,Covariates,Mortaility, Weighted Adj
table_dat <- X_df <- as.data.frame(cbind(X2)) %>%
  inner_join(Covariate_D[,c("SEQN","SDMVPSU","SDMVSTRA","WTMEC2YR")],by = "SEQN")%>%
  inner_join(Mortality_2011_D, by = "SEQN")%>%
  mutate(weight_adj = WTMEC2YR/mean(WTMEC2YR)) %>%
  mutate( yr5_mort =  as.integer(ifelse(permth_exm/12 <= 5 & mortstat == 1, 1,
                                        ifelse(permth_exm/12 < 5 & mortstat == 0, NA, 0))))%>%
filter(complete.cases(yr5_mort))
# #Train/Test Split
 split.train <- c(sample(which(X_df$yr5_mort ==0), 0.7*sum(X_df$yr5_mort ==0)),
                 sample(which(X_df$yr5_mort ==1),0.7*sum(X_df$yr5_mort == 1)))
 split.test <- setdiff(1:nrow(X_df), split.train)
  X_train <- X_df[split.train,]
  X_test <- X_df[split.test,]

#---------------------------------------------------------
#Predict mortality
predict.mortality <- function(scores_m,type = "Logistic",AIC = FALSE,train = X_train, test = X_test){
X.train <- inner_join(train,as.data.frame(scores_m), by = "SEQN")
X.test <- inner_join(test,as.data.frame(scores_m), by = "SEQN")
X <- inner_join(X_df,as.data.frame(scores_m), by = "SEQN" )
## Create svydesign() objects
X_svy.train <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight_adj, data = X.train, nest = TRUE, probs = NULL)
X_svy.test <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight_adj, data = X.test, nest = TRUE, probs = NULL)
X_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight_adj, data = X, nest = TRUE, probs = NULL)

#---------------------------------------------------------------
#AIC Selection
ind_vars  <- colnames(scores_m)[-1]
inc_vars  <-setdiff(names(X2), c("SEQN","intercept"))
exc_vars <- ind_vars
aic_vec  <- var_vec <- model_vec <- rep(NA, length(ind_vars))
## get the final model as the first model where AIC increases after removing a variable
backward_model <- paste0(c(ind_vars, inc_vars),collapse = "+")
if(AIC == TRUE){
  for(i in 1:length(ind_vars)){
    aic_ij  <- rep(NA,length(exc_vars))
    
    for(k in 1:length(exc_vars)){
      form    <- paste0(c(inc_vars, exc_vars[-k]), collapse="+")
      fit_tmp <- svyglm(as.formula(paste("yr5_mort ~", form)), design=X_svy.train,family=binomial())
      aic_ij[k] <- fit_tmp$aic
      rm(list=c("fit_tmp","form"))
    }
    
    k_cur         <- which(aic_ij == min(aic_ij))
    model_vec[i]  <- paste0(c(inc_vars, exc_vars[-k_cur]), collapse="+")
    exc_vars      <- exc_vars[-k_cur]
    aic_vec[i]    <- aic_ij[k_cur]
    rm(list=c("k_cur","aic_ij","k"))
  }
  backward_model <- model_vec[which(diff(aic_vec) > 0) + 1][1]
  if(is.na(backward_model))backward_model <- model_vec[1][1]
}
vars_select <<- strsplit(backward_model, "\\+")[[1]]

#---------------------------------------------
## estimating survey-weighted generalized linear models
if(type == "Logistic"){
mod.fit <- svyglm(as.formula(paste("yr5_mort ~", backward_model)), design=X_svy.train,  family=quasibinomial())
y_hat <- predict(mod.fit,X.test[,c(vars_select)],weights = X.test$weight_adj ,type = "response")
##Weighted ROC Curves, AUC
pred_ROCR = WeightedROC(y_hat, X.test$yr5_mort, weight = X.test$weight_adj)
}

## estimating survey-weighted boosting model.
if(type == "Boosting"){
  mod.fit <<- gbm(formula = as.formula(paste("yr5_mort ~", backward_model)), data = X_svy.train$variables,weights = X_train$weight_adj,
                  distribution = "adaboost",n.trees = 500,cv.folds = 10)
  y_hat <<- predict(mod.fit,newdata = X.test[,c(vars_select)],weights = X.test$weight_adj, type = "response") 
  pred_ROCR = WeightedROC(y_hat, X.test$yr5_mort, weight = X.test$weight_adj)
}

## estimating survey-weighted random forest model.
if(type == "RandomForest"){
  mod.fit <<- ranger(as.formula(paste("yr5_mort ~", backward_model)), data = X.train, case.weights = X.train$weight_adj)
  y_hat <<- predict(mod.fit, X.test[,c(vars_select)], case.weights = X.test$weight_adj)$predictions
  pred_ROCR = WeightedROC(y_hat, X.test$yr5_mort, weight = X.test$weight_adj)
}
AUC = WeightedAUC(pred_ROCR)
return(list(model_fit = mod.fit, yhat = y_hat, roc_ROCR = pred_ROCR,AUC = AUC))
}

#Classification Results:Logistic
narfd_logistic <- predict.mortality(scores_narfd,AIC = TRUE)
pfpca_logistic <- predict.mortality(scores_pfpca,AIC = TRUE)
gfpca_logistic <- predict.mortality(scores_gauss,AIC = TRUE)
nofpca_logistic <- predict.mortality(scores_m = data.frame(SEQN = seqn),AIC = FALSE)

#Classification Results:Boosting
narfd_gbm <- predict.mortality(scores_narfd, type = "Boosting")
pfpca_gbm <- predict.mortality(scores_pfpca, type = "Boosting")
gfpca_gbm <- predict.mortality(scores_gauss, type = "Boosting")
nofpca_gbm <- predict.mortality(scores_m = data.frame(SEQN = seqn), type = "Boosting",AIC = FALSE)

#Classification Results:Random Forest
narfd_rf <- predict.mortality(scores_narfd, type = "RandomForest")
pfpca_rf <- predict.mortality(scores_pfpca, type = "RandomForest")
gfpca_rf <- predict.mortality(scores_gauss, type = "RandomForest")
nofpca_rf <- predict.mortality(scores_m = data.frame(SEQN = seqn), type = "RandomForest",AIC = FALSE)

#Regression Fits: Logistic
narfd_logistic2 <- predict.mortality(scores_narfd,train = X_df, test = X_df,AIC = TRUE)
pfpca_logistic2 <- predict.mortality(scores_pfpca,train = X_df, test = X_df,AIC = TRUE)
gfpca_logistic2 <- predict.mortality(scores_gauss,train = X_df, test = X_df,AIC = TRUE)

rm(list = c('scores_gauss', 'scores_narfd' , 'scores_pfpca','predict.mortality'))
rm(list = c("X_svy","backward_model","X2","X_df","ind_vars","inc_vars","exc_vars"))
