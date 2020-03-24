###############################################################################
# If not already done...
###############################################################################

# Install relevant packages of 'ggtree' family
# First, BiocManager has to installed from https://bioconductor.org/install/
# Accordind to instructions there, use following commands:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

# install necessary packages from BiocManager
# For me, to have this working, removing  
# the previously installed Rccp package was necessary.
# Also, I agreed to update some of the packages when this was asked in the console 
BiocManager::install(c("ggtree", "treeio"))

###############################################################################
# Load required packages
###############################################################################

library(here)
library(ape)
library(ggtree)
library(ggstance)
library(Hmisc)
library(gdata) # messages on absent perl interpreter can be ignored
library(treeio)

# Set the current directory as working directory
set_here()


###############################################################################
# Load data
###############################################################################

#Load tree 
tree.unr <- read.tree("perl_its_PhyML_newick_tree.nhx")
#to check the tree struture use
#str(tree.unr) 


#Rooting and adjusting tip order
outgroup <- c("Ordynets_00158", "TU_124387",  
              "TU_124388", "ARAN_4160", "KAS_L_0103")
tree <-root(tree.unr, outgroup)


# Loading metadata to use for tree annotation
annot.full <- read.csv ("SUB_cystidia_measurements_20190808_v1.2.csv",  sep=",")
annot.perl1 <- annot.full[which(annot.full$Species_ID=="perlongisporum"), ]
annot.perl <- annot.perl1[c("Specimen_ID", "L_sp", "W_sp", "Q_sp")]



###############################################################################
# Visualize the tree with annotations as boxplots
###############################################################################

# Generate treedata object
p0<-ggtree(tree)


p1 <- p0 + 
  geom_tree(size=0.7) +
  geom_tiplab(size=2.5, align=F, hjust=-0.01)


                
p2<- facet_plot(p1+xlim_tree(0.25), panel='Cystidia length, mkm', data=annot.perl, 
                geom=geom_boxploth, mapping = aes(x=L_sp, group=label)) 

p3<- facet_plot(p2, panel='Cystidia width, mkm', data=annot.perl, 
                geom=geom_boxploth, mapping = aes(x=W_sp, group=label)) 

p4<- facet_plot(p3, panel='Cystidia length/width ratio', data=annot.perl, 
                geom=geom_boxploth, mapping = aes(x=Q_sp, group=label)) +
  theme_tree2()

print(p4)

# End