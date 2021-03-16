## In order to perform survey weighted PCA on the individual daily profiles,
## we first calculate weights for each "day" by dividing and individual's 4 year
## adjusted survey weight by the number of days of data they have.
## These daily weights are then used in the survey weighted pca calculation.
uid   <- unique(data_analysis$SEQN)
nid   <- length(uid)
n_uid <- vapply(uid, function(x) sum(Act_Analysis$SEQN == x), numeric(1))
wtmec4yr_adj_norm_rep <- SDMVPSU_rep <- SDMVSTRA_rep <- c()
for(i in 1:nid){
    inx                   <- which(data_analysis$SEQN == uid[i])
    wtmec4yr_adj_norm_rep <- c(wtmec4yr_adj_norm_rep, rep(data_analysis$wtmec4yr_adj_norm[inx]/n_uid[i],n_uid[i]) )
    SDMVPSU_rep           <- c(SDMVPSU_rep, rep(data_analysis$SDMVPSU[inx],n_uid[i]))
    SDMVSTRA_rep          <- c(SDMVSTRA_rep, rep(data_analysis$SDMVSTRA[inx],n_uid[i]))
}

data_svypca <- data.frame(Act,
                          "wtmec4yr_adj_norm_rep" = wtmec4yr_adj_norm_rep,
                          "SDMVPSU_rep" = SDMVPSU_rep,
                          "SDMVSTRA_rep" = SDMVSTRA_rep)
data_svypca <- svydesign(id= ~SDMVPSU_rep, strata = ~SDMVSTRA_rep, weights = ~wtmec4yr_adj_norm_rep,
                         data = data_svypca, nest = TRUE)
svypca_fit <- svyprcomp(formula=as.formula(paste0("~", paste0("MIN",1:1440, collapse="+") )),
                        design=data_svypca, scores=TRUE)




textscale <- 2
jpeg(file.path(figure_path, "functional_principal_components_vs_svyprcomp.jpeg"),height=1350,width=1350,quality=100)
par(mfrow=c(4,4))
for(k in 1:16){
    x <- svypca_fit$rotation[,k]
    y <- fpca_fit$efunctions[,k]
    matplot(1:1440,cbind(x,y), type="l", col=c("red","black"),xaxt='n',xlab="Time of Day", main=paste0("PC ", k),
            ylab="log(1+AC)", cex.main=textscale,cex.lab=textscale, cex.axis=textscale,
            lwd=1.5*textscale,cex=1.5*textscale)
    axis(1, at=c(1,6,12,18,23)*60 +1 , labels=c("01:00","06:00","12:00","18:00","23:00"),
         cex=textscale,cex.lab=textscale, cex.axis=textscale)



    tot_svy_inx <- min(which(cumsum(svypca_fit$sdev^2)/sum(svypca_fit$sdev^2) > 0.99))
    fpca_pct_var <- round(100*fpca_fit$evalues[k]/sum(fpca_fit$evalues),1)
    svy_pct_var <- round(100*svypca_fit$sdev[k]^2/sum(svypca_fit$sdev[1:tot_svy_inx]^2),1)

    legend("topright", c(paste0("svyprcomp: ", svy_pct_var, " % Explained"),
                        paste0("fPCA: ", fpca_pct_var, " % Explained")),
           col=c("red","black"),lty=1,lwd=textscale,bty='n',cex=textscale)

    rm(list=c("x","y","svy_pct_var","fpca_pct_var","tot_svy_inx"))
}
dev.off()



rm(list=c("uid","nid","i","data_svypca","k","textscale","n_uid",
          "svypca_fit","SDMVSTRA_rep","SDMVPSU_rep","wtmec4yr_adj_norm_rep","inx"))





