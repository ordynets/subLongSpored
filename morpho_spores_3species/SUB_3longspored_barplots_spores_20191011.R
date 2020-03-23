#Reset R's brain
rm(list=ls())

#setwd tells R where to look
setwd("G:/R/SUB_perlongisporum/morpho v20191011")

#use getwd to confirm that R is now looking here
getwd()

#Load spore table and adjust it
sub <- read.csv ("Spore_ranges_specimens_20191011_v2.csv",  sep=",")
sub <- sub[c("Species_ID", "Lmean", "Wmean", "Qmean")]



# Quantile function
quantiles_50 <- function(x) {
  r <- quantile(x, probs=c(0.00, 0.25, 0.5, 0.75, 1.00))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


library(Hmisc)
library(reshape2)
library(ggplot2)

###############################################################
## species with fusiform spores
###############################################################

# Subsetting data for particular plots, reordering, renaming columns
long.lev <-c("cochleum", "longisporum", "perlongisporum")
sub$Species_ID <- factor(sub$Species_ID, levels=long.lev) #reordering factor levels
colnames(sub)<-c("Species", "Spore length", "Spore width", "Spore length/\n width ratio")

# Melting 
sub.melt<-melt(sub)
colnames(sub.melt)<-c("Species", "Variable", "Value")


# Plotting
ggplot(sub.melt, aes(x = Species, y = Value, fill=Variable, width=0.9))+
  stat_summary(fun.data = quantiles_50, geom="boxplot", position=position_dodge(1))+ 
  facet_grid(.~Species, scales="free")+
  labs(title="Basidiospores", x="Species", y="Absolute size (µm) or length/width ratio")+
  theme(text = element_text(size=20), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"), 
        axis.title.x = element_text(size=17, face="bold"),
        axis.title.y = element_text(size=17, face="bold", vjust=2),
        axis.text.y=element_text(colour="black"))+
  scale_y_continuous(breaks=seq(0,25, by=5), minor_breaks=NULL)+
  theme(legend.position="bottom")

#end