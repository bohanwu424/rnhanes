#Clear Workspace
dev.off()
cat("\014")
rm(list=ls())
set.seed(18552)


library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(reshape2)
library(rnhanesdata)
library(refund)

load("pa_nhanes_5min.RData")

#-----------------------------------------------
# Actigraphy count data: raw totals, for every 5 minute intervals
Ytot.m = melt(Ytot, varnames = c("Subject","Time"), value.name = "Ytot") %>% 
  mutate(
  t = sapply(as.character(Time),function(s)strsplit(s," ")[[1]][1]) %>% as.numeric()
)
  
y_1 = Ytot.m %>% filter(Subject == 2) 

p <- ggplot(y_1, aes(x=t, y=Ytot)) +
  geom_line( color="#69b3a2") + 
  xlab("") +
  theme_ipsum() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) 

#-----------------------------------------------
#Gaussian FPCA
Y_log = log(Y + 1)
fit.gaus = fpca.face(Y_log,npc = 5)

#Plot Egenfuncitons
plot_F(fit.gaus$efunctions)

rnhanesdata::Mortality_2011_C %>% dim()
