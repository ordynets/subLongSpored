
###############################################################################
# Start
###############################################################################

#Reset R's brain
rm(list=ls())

#setwd tells R in which directory to work with data
setwd("G:/R/SUB_perlongisporum/ggtree")

#use getwd to confirm that R is now looking here
getwd()



###############################################################################
# Load required packages
###############################################################################

library(ggtree)
library(treeio)
library(ggstance)
library(Hmisc)
library(gdata) # messages on absent perl interpreter can be ignored



###############################################################################
# Load data
###############################################################################

#Load tree 
tree <- read.tree("ITS perl fasttree")
#to check the tree struture use
#str(tree) 



# Loading metadata to use for tree annotation
annot.full <- read.csv ("SUB_spore_measurements_20190723_v1.csv",  sep=",")



###############################################################################
# Visualize the tree with annotations as violin plots
###############################################################################

#creat a tree as an R object 

p1<-ggtree(tree, aes(x, y)) + 
  geom_tree() + 
  theme_tree() + 
  #geom_treescale()+ 
  geom_tiplab(size=2, align=TRUE) 
#ignore warning message: Duplicated aesthetics after name standardisation: size 
#for more parameters see ggtree vignette for tree visualization

#tree can be checked via command 
#print(p2)


p2<- facet_plot(p1+xlim_tree(0.2), panel='Spore length, mkm', data=annot.full, 
           geom=geom_violinh, mapping = aes(x=L_sp, group=label)) 

p3<- facet_plot(p2, panel='Spore width, mkm', data=annot.full, 
                geom=geom_violinh, mapping = aes(x=W_sp, group=label)) 

p4<- facet_plot(p3, panel='Spore length/width ratio', data=annot.full, 
                geom=geom_violinh, mapping = aes(x=Q_sp, group=label)) +
  theme_tree2()

print(p4)

