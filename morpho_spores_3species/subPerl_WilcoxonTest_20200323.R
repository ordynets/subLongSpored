# Load required packages
library(here)


# Set the current directory as working directory
set_here()


# Load spore table and adjust it
ranges.spo <- read.csv ("Spore_ranges_specimens_20191011_v2.csv",  sep=",")


# Get data for Subulocystidium perlongisporum only
ranges.spo.perl <- ranges.spo[which(ranges.spo$Species_ID=="perlongisporum"),]


# Mann-Whitney, = Wilcoxon tests
wil.sp.L<-wilcox.test(Lmean~perl_clade, data=ranges.spo.perl, lexact = FALSE)
wil.sp.L
wil.sp.W<-wilcox.test(Wmean~perl_clade, data=ranges.spo.perl,  exact = FALSE)
wil.sp.W
wil.sp.Q<-wilcox.test(Qmean~perl_clade, data=ranges.spo.perl,  exact = FALSE)
wil.sp.Q

# End
