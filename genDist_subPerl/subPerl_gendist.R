
###############################################################################
# Load required packages
###############################################################################
library(here)
library(ape)
library(plyr)
library(dplyr)
library(adegenet)
library(ggplot2)
library(forcats)


#Set the current directory as working directory
set_here()


###############################################################################
# Load DNA sequence alignment and analyse it 
###############################################################################
#ITS alignment
its <-read.dna("perl_ITS_27x_20190725_mafft_ITSx.phy")

# Plot alignment with legend, just for quick visual check
par(mar=c(3,5,3,3)) 
image(its, cex=0.5)


#Check naming and ordering of sequences in the alignment
labels(its)


#Assign sequences to clade 1 or 2 of Subulicystidium perlongisporum
#In this case, I can made it so easy because
#the sequences from two clades are not intermixed
grp<-c(rep("clade 1", 22), rep("clade 2", 5))


#Aplly function from "adegenet" to get within/between group distances
#distance method is set as in ape4s gendist
tempT <- pairDistPlot(its, grp, model="raw", pairwise.deletion = TRUE,
                     within=TRUE, sep="-", data=TRUE,
                     violin=TRUE, boxplot=TRUE, jitter=TRUE)


#Have a look onto output of pairDistPlot()
tempT$data
tempT$boxplot
tempT$violin


groupstat<-tempT$data %>%
  group_by(groups) %>%
  summarise(mindist=min(distance),
            maxdist=max(distance),
            mediandist=median(distance))

#Change ordering of levels in tempT$data
#and create a separate object for this dataframe
p <-tempT$data %>%
  mutate(groups = fct_relevel(groups, "clade 1-clade 1", "clade 2-clade 2", "clade 1-clade 2")) 


#Change naming of levels for more convenient plotting
levels(p$groups)[levels(p$groups)=="clade 1-clade 1"] <- "Within \n clade 1"
levels(p$groups)[levels(p$groups)=="clade 2-clade 2"] <- "Within \n clade 2"
levels(p$groups)[levels(p$groups)=="clade 1-clade 2"] <- "Clade 1 vs \n clade 2"


#Plot distances with ggplot
ggplot(p, aes(groups, distance)) + 
  geom_jitter(aes(colour=groups), 
              width = 0.25, size=6.0, shape=1, fill = NA, stroke = 2, alpha=0.5)+
  labs(x="Clades", y="Pairwise uncorrected genetic distances")+
  theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=17, face="bold", vjust=3.5),
        axis.text.x=element_text(colour="grey40", size=15),
        axis.text.y=element_text(colour="grey40", size=15))+ 
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.6,0.6,0.6,0.6), "cm"))

###############################################################################
# Write dataframe output from pairDistPlot to csv file
###############################################################################

write.table(tempT$data, file="SUB_perl_withinBetweenCladeDist_YYYYMMDD_v0.csv",sep=",",quote=F)

