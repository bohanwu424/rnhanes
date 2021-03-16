source(here::here("source", "NARFD_Utilities.R"))

# inputs:
#   long--data.frame with columns .index (time value), .id (curve label), 
#         .observed (observation of curve at that time value)
#   npc--how many FPCs/functional prototypes to estimate
#   nbasis--how many spline basis functions to use
#   type--NARFD of GFPCA
#   D--how many knots to use in spline basis
#   folds--integer, how many folds to use in cross-validation
#   iter--integer, maximum number of iterations to be used in alternating minimization
#   periodic--T if the spline basis used should be periodic, F otherwise
#   penalty.seq--numeric vector of penalty values to try
#   verbose--T to print out information about progress

# returns a list with elements
#   pred: data.frame with columns 
#         .index--time value
#         .id--curve label
#         .observed--count data
#         .pred--predicted mean (for GFPCA, after exponentiation)
#          prototype1, prototype2, ...--contributions to predicted mean from each prototype
#             (for GFPCA, prior to exponentiation)
#   scores: data.frame with columns prototype, .id and Value 
#           (Value is the estimated score for that prototype and curve)
#   prototypes: data.frame with columns prototype, .index and Value
#           (Value is the value of the prototype at the corresponding .index)
#   penalty: the value of the selected penalty
#   mean: for NARFD, NA, 
#     for GFPCA, a data.frame with columns .index and Value
#      
NARFD <- function(long, npc, nbasis, type, D, 
                    folds, iter, periodic, 
                    penalty.seq, verbose = F){
  subjs <- unique(long$.id)
  n.curves <- length(subjs)

  # get spline basis functions
  if (periodic) bs.func <- pbs else bs.func <- bs
  Theta <- bs.func(1:D, df=nbasis, intercept=T)
  # this is for fast access to the spline basis
  Theta.dt <- GetThetaDt(Theta, unique.points=unique(long$.index))
  # get initial values
  inits <- GetInits(type=type, long=long, npc=npc, Theta=Theta)
  # set up folds
  fold.df <- data.table(Curve=subjs, 
                        Which.fold=sample(rep(1:folds, ceiling(n.curves / folds))[1:n.curves], 
                                          size=n.curves, 
                                          replace=F))
  long <- as.data.table(long)
  
  all.errors <- list()
  for (f in 1:folds){
    if (verbose) {cat("Doing fold", f, "\n")}
    for (i in 1:length(penalty.seq)){
      p <- penalty.seq[i]
      if (verbose) {cat("penalty", p, "\n")}
      r <- CrossValidate(fold.df=fold.df, f=f, p=p, long=long, 
                         npc=npc, inits=inits, Theta=Theta, 
                         type=type, iter=iter, Theta.dt=Theta.dt)
      this.error <- r$error
      if (i > 1){
        fold.errors <- rbind(fold.errors, this.error)
      } else {
        fold.errors <- this.error
      }
    }
    all.errors[[f]] <- fold.errors
  }
  errors <- bind_rows(all.errors)
  
  error.summ <- group_by(errors, Penalty) %>% summarize(Mean.MKL=mean(mkl)) %>% arrange(Mean.MKL)
  best.penalty <- filter(error.summ, Mean.MKL==min(Mean.MKL))$Penalty
  # once the best penalty is fixed, redo with all the data and 
  # the selected value
  res.l <- AlternatingFit(dat=long, npc=npc, 
                          inits=inits,
                          Theta=Theta, 
                          type=type, 
                          iter=iter,
                          lambda=best.penalty, 
                          decompose=T, verbose = verbose)
  pred <- res.l$pred %>% as.data.frame()
  pcs <- ProcessPCs(Theta %*% t(res.l$pc.coefs))
  scores <- ProcessScores(res.l$scores)
  l <- ReorderAndScale(pred, scores, pcs, type=type)
  if (type=="GFPCA"){
    mean <- data.frame(.index=1:nrow(Theta), 
                       Value=Theta %*% t(res.l$mean.coefs))
  } else {
    mean <- NA
  }
  # do PCs, do scores, do mean, return pred, return pen.summ and pen
  colnames(l$pred) <- gsub("FPC", "prototype", colnames(l$pred))
  l$scores <- l$scores %>%
    rename(prototype = FPC) %>%
    mutate(prototype = gsub("FPC", "prototype", prototype))
  l$pcs <- l$pcs %>%
    rename(prototype = FPC) %>%
    mutate(prototype = gsub("FPC", "prototype", prototype))
  list(pred=l$pred, scores=l$scores, prototypes=l$pcs, 
       penalty=best.penalty, 
       mean=mean)
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
