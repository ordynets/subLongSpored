

#Reset R's brain
rm(list=ls())

#setwd tells R where to look
setwd("G:/R/SUB_perlongisporum/morpho v20191011")

#use getwd to confirm that R is now looking here
getwd()

#Load spore table and adjust it
ranges.spo <- read.csv ("Spore_ranges_specimens_20191011_v2.csv",  sep=",")


ranges.spo.perl <- ranges.spo[which(ranges.spo$Species_ID=="perlongisporum"),]


# Mann-Whitney, = Wilcoxon tests
wil.sp.L<-wilcox.test(Lmean~perl_clade, data=ranges.spo.perl, lexact = FALSE)
wil.sp.L
wil.sp.W<-wilcox.test(Wmean~perl_clade, data=ranges.spo.perl,  exact = FALSE)
wil.sp.W
wil.sp.Q<-wilcox.test(Qmean~perl_clade, data=ranges.spo.perl,  exact = FALSE)
wil.sp.Q
