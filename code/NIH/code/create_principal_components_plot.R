sd_pcs <- apply(fpca_fit$scores,2,sd)[1:8]
inx_plot <- c(1:1440)[seq(1,1440,by=20)]
textscale <- 2
jpeg(file.path(figure_path, "functional_principal_components.jpeg"),height=900,width=1350,quality=100)
par(mfrow=c(2,3),las=1,mar=c(5.1,6.1,2.1,2.1))
for(i in 1:6){
    low_tmp <- fpca_fit$mu + -2*sd_pcs[i]*fpca_fit$efunctions[,i]
    high_tmp <- fpca_fit$mu + 2*sd_pcs[i]*fpca_fit$efunctions[,i]
    matplot(inx_plot,cbind(low_tmp,fpca_fit$mu,high_tmp)[inx_plot,], type=c('p',"l","p"), pch=c("-",NA,"+"),lty=c(1,1,1),
            col=c("grey","black","grey"),xaxt='n',xlab="Time of Day",ylim=c(-0.5,6),
            main=paste("PC ", i, " Percentage of variabiity ", round(100*fpca_fit$evalues[i]/sum(fpca_fit$evalues),1), "%", sep=""),
            ylab="log(1+AC)", cex.main=textscale,cex.lab=textscale, cex.axis=textscale,
            lwd=1.5*textscale,cex=1.5*textscale)
    axis(1, at=c(1,6,12,18,23)*60 +1 , labels=c("01:00","06:00","12:00","18:00","23:00"),
         cex=textscale,cex.lab=textscale, cex.axis=textscale)
}
dev.off()
rm(list=c("sd_pcs","inx_plot","high_tmp","low_tmp","i","textscale"))






