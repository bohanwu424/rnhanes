library(dplyr)
plot_F <- function(F){
  F_d = melt(F)  %>% dplyr::rename(.index = Var1,prototype = Var2, Value = value)
  p = F_d%>% 
    ggplot(aes(x = .index, y = Value)) + 
    geom_line() + 
    facet_grid(~prototype) + 
    xlab("Time index")
  return(p)
}
