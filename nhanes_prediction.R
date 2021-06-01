rm(list = ls())
{
options(
  ggplot2.continuous.colour = "viridis",
  ggplot2.continuous.fill = "viridis"
)

scale_colour_discrete = scale_colour_viridis_d
scale_fill_discrete = scale_fill_viridis_d

theme_set(theme_minimal() + theme(legend.position = "right"))
}
#Previously in the code "nhanes_feature_extraction", We fitted GFPCA, NARFD, and PFPCA
#to extract dominant functional features from Nhanes accelerometry data and visualize the results. 

#In this code,we construct prediction models using logistic regression, random forest and AdaBoost,
#combined with scores from each one of GFPCA,NARFD, and PFPCA model fits. 
#Before applying the prediction methods, we use regresion-based AIC selection to choose the set of predictors. 
#In total, we have nine prediction models to compare.
#The AIC selection, prediction modeling, and results are all survey-weighted.  
#Finally, the code provides useful visuals and intepretations for the prediction results.   

source(here::here("source", "Helper Functions.R"))
load(here::here("data", "NARFD_results.RData"))
load(here::here("data", "PFPCA_results.RData"))
load(here::here("data", "GFPCA_results.RData"))

#--------------------------------------------
## Import data

load(here::here("data", "pa_nhanes_5min.RData"))

rownames(Y) <- seq(nrow(Y)); colnames(Y) <- seq(ncol(Y))

#----------------------------------------
#prediction
source(here::here("source", "predict_mortality.R"))

#----------------------------------------------
#AUC Table for Logistic, Random Forest, AdaBoost with GFPCA, PFPCA, and NARFD scores
AUC <- data.frame(model = rep(c("Logistic", "Random Forest","AdaBoost"),each = 4),
                  submodel = rep(c("GFPCA","PFPCA","NARFD","Baseline"),3), 
                  AUC = c(round(gfpca_logistic$AUC,3),round(pfpca_logistic$AUC,3), round(narfd_logistic$AUC,3),round(nofpca_logistic$AUC,3),
                          round(gfpca_rf$AUC,3),round(pfpca_rf$AUC,3), round(narfd_rf$AUC,3),round(nofpca_rf$AUC,3),
                          round(gfpca_gbm$AUC,3),round(pfpca_gbm$AUC,3), round(narfd_gbm$AUC,3),round(nofpca_gbm$AUC,3)))

AUC <- AUC %>% arrange(desc(AUC)) %>% 
  recast(formula = model ~ submodel)

table_auc <- kable(AUC,format="latex",booktabs=TRUE,longtable =TRUE, caption="AUC scores for mortality prediction") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
write(table_auc,file=file.path("Tables","table_auc.tex") )

#------------------------------------------------------
#Plot ROC Curve for all
preds = list(p1 = narfd_gbm$roc_ROCR, p2 = pfpca_gbm$roc_ROCR, p3 = gfpca_gbm$roc_ROCR,
             p4 = narfd_logistic$roc_ROCR, p5 = pfpca_logistic$roc_ROCR, p6 = gfpca_logistic$roc_ROCR, 
             p7 = narfd_rf$roc_ROCR, p8 = narfd_rf$roc_ROCR, p9 = narfd_rf$roc_ROCR)

jpeg(file.path("plots","ROC comparison.jpeg"),height=950,width=1400,quality = 100)

ggplot() + 
  geom_line(aes(preds[[1]][, c("FPR")], preds[[1]][, c("TPR")], color = "NARFD", linetype = "AdaBoost"),size = 3)+
  geom_line(aes(preds[[2]][, c("FPR")], preds[[2]][, c("TPR")],color = "Poisson FPCA",linetype = "AdaBoost"),size = 3)+
  geom_line(aes(preds[[3]][, c("FPR")], preds[[3]][, c("TPR")], color = "Gaussian FPCA",linetype = "AdaBoost"),size = 3)+
  geom_line(aes(preds[[4]][, c("FPR")], preds[[4]][, c("TPR")], color = "NARFD",linetype = "logistic"),size = 3)+
  geom_line(aes(preds[[5]][, c("FPR")], preds[[5]][, c("TPR")],color = "Poisson FPCA",linetype = "logistic"),size = 3)+
  geom_line(aes(preds[[6]][, c("FPR")], preds[[6]][, c("TPR")], color = "Gaussian FPCA",linetype = "logistic"),size = 3)+
  geom_line(aes(preds[[7]][, c("FPR")], preds[[7]][, c("TPR")], color = "NARFD",linetype = "random forest"),size = 3)+
  geom_line(aes(preds[[8]][, c("FPR")], preds[[8]][, c("TPR")],color = "Poisson FPCA",linetype = "random forest"),size = 3)+
  geom_line(aes(preds[[9]][, c("FPR")], preds[[9]][, c("TPR")], color = "Gaussian FPCA",linetype = "random forest"),size = 3)+
  scale_color_manual("", values=c("NARFD" = "blue", "Gaussian FPCA" = "green","Poisson FPCA" = "red"))+
  scale_linetype_manual("", values=c("random forest" = "dotted", "logistic" = "twodash","AdaBoost" = "solid"))+
  geom_abline(aes(slope = 1,intercept = 0))+
  labs(x = "False Positive Rate", y = "True Positive Rate")+
  ggtitle("ROC Curves for Mortality Prediction")+
  theme(text = element_text(size = 40))

dev.off()
#-----------------------------------------------------------
#Create Variable Importance Plots using XGBoost
{
varim.p <- summary(pfpca_gbm$model_fit)
varim.p$var = gsub("Education= HS","HS",varim.p$var)
varim.p$var = gsub("Above_HS","Education> HS",varim.p$var )
varim.p$var = gsub("Chol_HDL","`HDL Cholesterol`", varim.p$var)
varim.p$var = gsub("Chol_tot","`Total Cholesterol`", varim.p$var)
varim.p$var = gsub("Blood_Pressure","`Systolic Blood Pressure`", varim.p$var)
varim.p$var = gsub("n_wkday","`# of Weekdays Recorded`" , varim.p$var)
varim.p$var = gsub("n_wkend","`# of Weekend Days Recorded`", varim.p$var)
}

jpeg(file.path("plots","Var importance(AdaBoost + PFPCA).jpeg"),height=900,width=1350,quality = 100)
ggplot(varim.p, aes(x =reorder(var, rel.inf), y = rel.inf))+
  geom_col()+
  labs(x = "Variable",  y = "Rel importance")+
  ggtitle("Variable importance (AdaBoost + PFPCA)")+
  theme(text = element_text(size = 27))+
  labs(y = "Relative importance", x = "Variable")+
  coord_flip()
dev.off()
#------------------------------------------------
##Create Table 1
#source(here::here("source", "create_table1.R"))

#------------------------------------------------
##Create coefficient table for logistic regression
{
  coeff <- coefficients(pfpca_logistic2$model_fit)
  names(coeff) = gsub("Education= HS","HS",names(coeff))
  names(coeff) = gsub("Above_HS","Education> HS",names(coeff) )
  names(coeff) = gsub("Chol_HDL","`HDL Cholesterol`", names(coeff))
  names(coeff) = gsub("Chol_tot","`Total Cholesterol`", names(coeff))
  names(coeff) = gsub("Blood_Pressure","`Systolic Blood Pressure`", names(coeff))
  names(coeff) = gsub("n_wkday","`# of Weekdays Recorded`" , names(coeff))
  names(coeff) = gsub("n_wkend","`# of Weekend Days Recorded`", names(coeff))
  coeff_CI <- round(cbind(coeff,confint(pfpca_logistic2$model_fit)),3)
  coeff_CI <- data.frame(Variable = rownames(coeff_CI),
                         `Coefficient (95% Confidence Interval)` = apply(coeff_CI,1,function(x) paste0(x[1]," (",x[2],", ",x[3],")")))
  confint.svyglm
}
coeff_CI <- as.matrix(coeff_CI) ; rownames(coeff_CI) = c()
table_coef1 <- kable(coeff_CI,format="latex",booktabs=TRUE,longtable =FALSE, caption="Coefficents table for logistic + GFPCA") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
write(table_coef1,file=file.path("Tables","table_coef3.tex") )

rm(list = c("scores_m","samples","roc_ROCR","AUC"))
library(survey)
summary(pfpca_logistic2$model_fit)
