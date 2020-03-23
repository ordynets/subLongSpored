#Reset R's brain
rm(list=ls())

#setwd tells R where to look
setwd("G:/R/SUB_perlongisporum/morpho v20191011")

#use getwd to confirm that R is now looking here
getwd()

#Load spore table and adjust it
sub <- read.csv ("SUB_cystidia_measurements_20191011_v1.csv",  sep=",")



library(Hmisc)
library(dplyr)
library(tolerance)


# S. perlongisporum
# Subset relevant measuremnents
# Then calculate mean per specimen
perlongisporum <-sub[sub$Species_ID %in% "perlongisporum", ] 
perlongisporum.means <- perlongisporum %>% 
  group_by(Specimen_ID) %>% 
  summarise(Lm=mean(L_sp), Wm=mean(W_sp), Qm=mean(Q_sp))
perlongisporum.interv.L <-normtol.int(perlongisporum.means$Lm, alpha = 0.1, P = 0.9, side = 2, method = "HE")
perlongisporum.interv.L
perlongisporum.interv.W <-normtol.int(perlongisporum.means$Wm, alpha = 0.1, P = 0.9, side = 2, method = "HE")
perlongisporum.interv.W
perlongisporum.interv.Q <-normtol.int(perlongisporum.means$Qm, alpha = 0.1, P = 0.9, side = 2, method = "HE")
perlongisporum.interv.Q

# End of the script
