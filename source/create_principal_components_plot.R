plot_pc <- function(Y,scores_m, F_m,mu = NA, textscale = 2){
  sd_pcs <- apply(scores_m,2,sd)[1:8]
  T = nrow(F_m)
  inx_plot <- c(1:T)[seq(1,T,by=20)]
  jpeg(file.path("plots", "functional_principal_components.jpeg"),height=900,width=1350,quality=100)
  par(mfrow=c(2,3),las=1,mar=c(5.1,6.1,2.1,2.1))
  for(i in 1:6){
    if(is.na(mu)){
      mu = colMeans(Y)
    }
    low_tmp <-  mu -2*sd_pcs[i]*F_m[,i]
    high_tmp <- mu + 2*sd_pcs[i]*F_m[,i]
    matplot(inx_plot,cbind(low_tmp,mu,high_tmp)[inx_plot,], type=c('p',"l","p"), pch=c("-",NA,"+"),lty=c(1,1,1),
            col=c("darkgrey","black","darkgrey"),xaxt='n',xlab="Time of Day",
            main=paste("PC ", i, sep=""),
            ylab="AC", cex.main=textscale,cex.lab=textscale, cex.axis=textscale,
            lwd=1.5*textscale,cex=1.5*textscale)
    axis(1, at=c(1,6,12,18,23)*12 +1 , labels=c("01:00","06:00","12:00","18:00","23:00"),
         cex=textscale,cex.lab=textscale, cex.axis=textscale)
  }
  dev.off()
  rm(list=c("sd_pcs","inx_plot","high_tmp","low_tmp","i","textscale"))
}

plot_pc(Y,scores_m = scores_m, F_m = fpcs_m,textscale = 3)

rm(list = c("plot_pc"))



