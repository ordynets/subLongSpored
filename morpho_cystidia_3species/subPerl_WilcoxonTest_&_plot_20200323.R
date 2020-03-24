# Load required packages
library(here)
library(ggplot2)
library(ggstance)


#Load the input table and adjust it
ranges.cys <- read.csv ("Cystidia_ranges_specimens_20191011_v2.csv",  sep=",")


# Get the data for the single species, Subulicystidium perlongisporum
ranges.cys.perl <- ranges.cys[which(ranges.cys$Species_ID=="perlongisporum"),]


# Mann-Whitney, = Wilcoxon tests
# legitimate test in this case
wil.cy.L<-wilcox.test(Lmean~perl_clade, data=ranges.cys.perl, exact = FALSE)
wil.cy.L
wil.cy.W<-wilcox.test(Wmean~perl_clade, data=ranges.cys.perl,  exact = FALSE)
wil.cy.W
wil.cy.Q<-wilcox.test(Qmean~perl_clade, data=ranges.cys.perl,  exact = FALSE)
wil.cy.Q



# Plot boxplot with dotplots
# Outliers points are marked with a dot in the center via the function geom_boxplot()  

bp <- ggplot(ranges.cys.perl, aes(x=factor(perl_clade, level=c("clade 1", "clade 2")), 
                                  y=Lmean)) + 
  geom_point(size=5, shape=21, fill = NA, aes(x=factor(perl_clade, level=c("clade 1", "clade 2")), 
                                              y=Lmean, color=perl_clade))+
  geom_boxplot(fill = NA, outlier.size = 1, outlier.shape=19, aes(x=factor(perl_clade, level=c("clade 1", "clade 2")), 
                              y=Lmean, color=perl_clade))+
  labs(y="Mean length of cystidia, mkm", x="Clade of Subulicystidium perlongisporum")+
  theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title.x = element_text(size=14, face="bold", vjust=0.5),
        axis.title.y = element_text(size=14, face="bold", vjust=2.2),
        axis.text.x=element_text(colour="grey40", size=15),
        axis.text.y=element_text(colour="grey40", size=15))+
  theme(legend.position = "none")

bp

# End