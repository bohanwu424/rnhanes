


xlabs <- c("Age","ASTP","Cigarette","Gender","Alcohol","CHF","MVPA","Cancer","Surrogate: si6",
           "BMI","Diabetes","Education","Surrogate: si5","Surrogate: si1","Surrogate: mi2",
           "Mobility Problem","Stroke","CHD","Wear Time","SATP","TLAC","Sedentary Time","Surrogate: mi1", "Ethnicity")

actlab_inx <- which(xlabs %in% c("TLAC","SATP","ASTP","MVPA",
                                 "Sedentary Time", "Wear Time",
                                 "Surrogate: mi2","Surrogate: si1",
                                 "Surrogate: si6","Surrogate: si5","Surrogate: mi1"))
xlab_cols <- rep("black",length(xlabs)); xlab_cols[actlab_inx] <- "firebrick"
png("../AUCplot.png",height=1000,width=1600)
xinds <- 1:nrow(auc_mat_full)
par(mar=c(7,4.1,4.1,4.1),cex=2,las=1)
plot(xinds, auc_mat_full$AUC,xaxt="n",ylab="AUC",xlab="",main="Model Selection Criteria",type='b',pch=NA)

inx_auc <- which.min(diff(auc_mat_full[,2])>0)
inx_aic <- which.min(diff(auc_mat_full[,3])<0)
inx_epic <- which.min(diff(auc_mat_full[,4])<0)

cols <- rep("black",length(xinds))
cols[inx_auc] <- "firebrick"
text(xinds, auc_mat_full$AUC,labels=c("AUC"),cex=0.75,col=cols)
axis(1, at=xinds, labels=FALSE)
text(x=xinds, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
     labels=xlabs, srt=45, adj=1, xpd=TRUE,col=xlab_cols)
# axis(1,at=xinds,labels=xlabs,las=2)
par(new=TRUE)

ylims2 <- range(as.vector(unlist(auc_mat_full[,c("AIC","EPIC")])))
matplot(xinds,auc_mat_full[,c("AIC","EPIC")],axes=FALSE,xlab="",ylab="",main="",col=c("black","black"),type='b',lty=1,pch=c(NA,NA))
axis(side=4,at=pretty(ylims2),las=1)
mtext("AIC/EPIC", side=4,line=3,cex=2,las=0)

cols <- rep("black",length(xinds))
cols[inx_aic] <- "firebrick"
text(xinds, auc_mat_full$AIC,labels=c("AIC"),cex=0.75,col=cols)

cols <- rep("black",length(xinds))
cols[inx_epic] <- "firebrick"
text(xinds, auc_mat_full$EPIC,labels=c("EPIC"),cex=0.75,col=cols)

dev.off()
