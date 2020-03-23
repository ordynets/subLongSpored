

#Reset R's brain
rm(list=ls())

#setwd tells R where to look
setwd("G:/R/SUB_perlongisporum/morpho v20191011")

#use getwd to confirm that R is now looking here
getwd()

#Load cysre table and adjust it
ranges.cys <- read.csv ("Cystidia_ranges_specimens_20191011_v2.csv",  sep=",")


ranges.cys.perl <- ranges.cys[which(ranges.cys$Species_ID=="perlongisporum"),]


# Mann-Whitney, = Wilcoxon tests
# legitimate test in this case
wil.cy.L<-wilcox.test(Lmean~perl_clade, data=ranges.cys.perl, exact = FALSE)
wil.cy.L
wil.cy.W<-wilcox.test(Wmean~perl_clade, data=ranges.cys.perl,  exact = FALSE)
wil.cy.W
wil.cy.Q<-wilcox.test(Qmean~perl_clade, data=ranges.cys.perl,  exact = FALSE)
wil.cy.Q



#Plot boxplot with dotplots
library(ggplot2)
library(ggstance)







+ 
  theme_classic()

bp




 + 
  theme_classic()

bp



