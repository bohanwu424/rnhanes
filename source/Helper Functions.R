#devtools::install_github('julia-wrobel/registr')
#devtools::install_github('linxihui/NNLM')
#2devtools::install_github("andrew-leroux/rnhanesdata")

pckgs <-c("tableone","knitr","kableExtra", "devtools","reshape2","tidyverse",
          "survey", "mgcv","refund","rnhanesdata","gbm",
          "ranger","WeightedROC","registr")
sapply(pckgs, function(x) if(!require(x,character.only=TRUE,quietly=TRUE)) {
  install.packages(x)
  require(x, character.only=TRUE)
})
rm(list=c("pckgs"))

plot_F <- function(F){
  F_d = melt(F)  %>% dplyr::rename(.index = Var1,prototype = Var2, Value = value)
  p = F_d%>% 
    ggplot(aes(x = .index, y = Value)) + 
    geom_line() + 
    facet_grid(~prototype) + 
    xlab("Time index")
  return(p)
}

# converts wide data frame with observations in columns
# and converts to long data frame with columns .id, .observed and .index
# .id is character-valued, starting with I
# .index is the time index, starting with 1 and going to number of rows
# in wide
MakeLong <- function(wide){
  D <- nrow(wide)
  
  long <- as.data.frame(wide) %>% 
    mutate(.index = 1:D) %>% 
    gather(.id, .observed, -.index) %>%
    mutate(.id = paste0("I", .id)) %>% 
    as_tibble()
  
  long
}
MakeLong2 <- function(wide){
  D <- nrow(wide)
  
  long <- as.data.frame(wide) %>% 
    mutate(index = 1:D) %>% 
    gather(id, value, -index) %>%
    mutate(id = paste0("I", id)) %>% 
    as_tibble()
  
  long
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