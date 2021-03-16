library(rnhanesdata)
library(refund)
library(glmnet)
library(ROCR)
library(caret)
load("pa_nhanes_5min2.RData")


fpca_fit = fpca.face(log(1+as.matrix(Y)), knots = 50)

X = as.data.frame(X)

X[,paste0("FPC",1:6)] = fpca_fit$scores[,1:6]

X[,"SEQN"] = seqn

X = left_join(X, Mortality_2011_D, by = "SEQN")

X[,"yr5_mort"] <- as.integer(ifelse(X$permth_exm/12 <= 5 & X$mortstat == 1, 1,
                                       ifelse(X$permth_exm/12 < 5 & X$mortstat == 0, NA, 0)))

X = X[complete.cases(X[,"yr5_mort"]),]

test_fit = glm(yr5_mort ~., data =  X[,-(33:46)], family = "binomial")

pred_ROCR = prediction(test_fit$fitted.values, test_fit$y)

roc_ROCR = performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR)
abline(a=0, b=1)
auc_ROCR = performance(pred_ROCR, measure = "auc")
auc_ROCR = auc_ROCR@y.values[[1]]

confusionMatrix(as.factor(1*test_fit$fitted.values>0.5), as.factor(test_fit$y==1))
