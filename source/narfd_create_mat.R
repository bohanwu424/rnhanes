#Create Score/FPC matrices
make_mat <- function(score_m = pfpca_results$scores){
  scores <- cbind(seqn,score_m)
  colnames(scores) <- c("SEQN",paste0("score",seq(ncol(scores)-1)))
  scores
}
##PFPCA
scores_pfpca <- make_mat(pfpca_results$scores)
fpcs_pfpca = pfpca_results$efunctions
colnames(fpcs_pfpca) <- paste0("FPC", seq(ncol(fpcs_pfpca)))

##NARFD
scores_narfd <- reshape2::dcast(narfd_results$scores, .id ~ prototype,value.var = "Value") %>% 
  column_to_rownames(".id")%>%
  make_mat()
fpcs_narfd <- reshape2::dcast(narfd_results$prototypes, .index ~ prototype,value.var = "Value") %>% column_to_rownames(".index")
 colnames(fpcs_narfd) <- paste0("FPC", seq(ncol(fpcs_narfd)))

##Gaussian FPCA
scores_gauss <- make_mat(fpca_fit$scores[,1:6])
fpcs_gfpca = fpca_fit$efunctions
colnames(fpcs_gfpca) <- paste0("FPC", seq(ncol(fpcs_gfpca)))

# scores = list(gauss = scores_gauss, pfpca = scores_pfpca, narfd = scores_narfd)
# fpcs = list(gauss = fpcs_gfpca, pfpca = fpcs_pfpca, narfd = fpcs_narfd)
# mu = list(gauss = fpca_fit$mu, pfpca = pfpca_results$mu)
# save(scores,fpcs,mu ,file = "scores and fpcs.RData")
rm( list = c("make_mat","fpcs_pfpca","fpcs_narfd","fpcs_gfpca"))

