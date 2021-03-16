## get contributions for subjects 21039, 21051
## from the functional predictor
# subjs <- c(21039, 21051)
subjs <- c(21009,  21068)
lp    <- predict(fit_SoFR, newdata = subset(data_analysis, SEQN %in% subjs), type='iterms')
coef_tmp <- coef(fit_SoFR,n=1440)

lps <- sprintf("%5.2f",round(lp[,13],2))

ylims_act <- c(-0.75,7.5)
ylims_coef <- c(-3.25,3.25)

jpeg(file.path(figure_path,"estimated_functional_coefficient.jpeg"),height=1400,width=2800,quality=100)
layout(matrix(c(1,1,1,1,
                2,3,
                4,5),nrow=2,byrow=FALSE))
par(las=1,cex=2.5, mar=c(5.1, 3.1, 4.1, 3.1))
matplot(coef_tmp$X_cn.argvals,
        cbind(coef_tmp$value,coef_tmp$value + 1.96*coef_tmp$se, coef_tmp$value - 1.96*coef_tmp$se ),
        xlab="Time of the Day", ylab="", main =  expression(paste("(A) Functional Coefficient ",hat(gamma)(s))),
        xaxt='n',type='l', lty=c(1,2,2),col='black',lwd=2)
axis(1, at=(c(1,6,12,18,23)*60 +1)/1440 , labels=c("01:00","06:00","12:00","18:00","23:00"))
abline(h=0,col="slategrey",lty=2,lwd=6)

i1_inx <- which(Act_Analysis$SEQN %in% subjs[1])
i2_inx <- which(Act_Analysis$SEQN %in% subjs[2])

act_sm_tmp1 <- t(fpca_fit$Yhat[i1_inx,])
act_sm_tmp2 <- t(fpca_fit$Yhat[i2_inx,])
act_avg    <- colMeans(avg_profiles)
act_avg_sm <- cbind(rowMeans(act_sm_tmp1), rowMeans(act_sm_tmp2))
act_avg_sm_cn <- act_avg_sm - (act_avg %*% t(c(1,1)))


matplot(1:1440, act_sm_tmp1, pch=16,col=rgb(0,0,0,0.4),cex=0.5,lty=1,
        xlab='Time of the Day',ylab="log(1+AC)",main=bquote("(B) SEQN " ~ .(subjs[1]) ~ ": " ~ frac(1, J[i]) ~ sum(tilde(X)[ij](s))),xaxt='n',
        ylim=ylims_act,type='l')
axis(1, at=c(1,12,23)*60 +1, labels=c("01:00","12:00","23:00"))
lines(1:1440, act_avg_sm[,1],lwd=3.5)
lines(1:1440, act_avg, col="firebrick2",lwd=4.5,lty=2.5)

matplot(1:1440, act_sm_tmp2, pch=16,col=rgb(0,0,0,0.4),cex=0.5,lty=1,
        xlab='Time of the Day',ylab="log(1+AC)",main=bquote("(C) SEQN " ~ .(subjs[2]) ~ ": " ~ frac(1, J[i]) ~ sum(tilde(X)[ij](s))),xaxt='n',
        ylim=ylims_act,type='l')
axis(1, at=c(1,12,23)*60 +1, labels=c("01:00","12:00","23:00"))
lines(1:1440, act_avg_sm[,2],lwd=3.5)
lines(1:1440, act_avg, col="firebrick2",lwd=4.5,lty=2.5)


xx <- 1:1440
plot(1:1440, act_avg_sm_cn[,1]*coef_tmp$value, pch=16,col=rgb(0,0,0,0.5),cex=0.5,type='l',
     xlab='Time of the Day',ylab="",main=bquote("(D) " ~ .(subjs[1]) ~ ": (" ~ frac(1, J[i]) ~ sum(tilde(X)[ij](s))-bar(X) ~ ")" ~  hat(gamma)(s) == .(lps[1])),xaxt='n',
     ylim=ylims_coef)
axis(1, at=c(1,12,23)*60 +1, labels=c("01:00","12:00","23:00"))

polygon(x=c(xx,rev(xx)),
        y=c(act_avg_sm_cn[,1]*coef_tmp$value, rep(0,length(xx))),
        col="lightgray")


plot(1:1440, act_avg_sm_cn[,2]*coef_tmp$value, pch=16,col=rgb(0,0,0,0.5),cex=0.5,type='l',
     xlab='Time of the Day',ylab="",main=bquote("(E) " ~ .(subjs[2]) ~ ": (" ~ frac(1, J[i]) ~ sum(tilde(X)[ij](s))-bar(X) ~ ")" ~  hat(gamma)(s) == .(lps[2])),xaxt='n',
     ylim=ylims_coef)
axis(1, at=c(1,12,23)*60 +1, labels=c("01:00","12:00","23:00"))

polygon(x=c(xx,rev(xx)),
        y=c(act_avg_sm_cn[,2]*coef_tmp$value, rep(0,length(xx))),
        col="lightgray")


dev.off()

par(mar=c(5.1, 4.1, 4.1, 2.1))

rm(list=c("subjs","lp","lps", "coef_tmp","act_sm_tmp1",
          "act_sm_tmp2","act_avg_sm","xx","ylims_act","ylims_coef",
          "act_avg_sm_cn","i1_inx","i2_inx"))
