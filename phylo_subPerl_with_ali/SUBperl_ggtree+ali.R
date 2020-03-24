###############################################################################
# Install if not already done...
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
library(ape)
library(gdata)
library(Hmisc)
library(treeio)
library(ggtree)
library(ggrepel)

# Set the current directory as working directory
set_here()

# Read in data 
annot.full <- read.csv ("SUB_data_processing_20190731.csv",  sep="\t")

###############################################################################
# Load and manage tree tip labels
###############################################################################


#Load tree 
tree.unr <- read.newick("perl_its_PhyML_newick_tree.nhx")
#to check the tree struture use
str(tree.unr)


#Rooting and adjusting tip order
outgroup <- c("Ordynets_00158", "TU_124387",  
              "TU_124388", "ARAN_4160", "KAS_L_0103")
tree <-root(tree.unr, outgroup)


# First, sample from big annotation table the rows that correspond to the tips 
# in the given tree
ml.tip.lab<-tree$tip.label
annot.filter.ml <- annot.full[which(annot.full$Specimen_ID %in% ml.tip.lab), ]
annot.filter.ml <- drop.levels(annot.filter.ml) # drop unused, "ghost" levels


# Ordering species labels by the order of sequences in alignment/distance matrix
annot.match.ml <- annot.filter.ml[match(ml.tip.lab, annot.filter.ml$Specimen_ID),]

# Check that labels in tree and in table annot.match are 
# the same and in the same order
identical(tree$tip.label, as.character(annot.match.ml$Specimen_ID))



#Replace original tip labels in the tree by the ones needed for publication
tree$tip.label<-as.character(annot.match.ml$ITS_tip_labels)



###############################################################################
# # Load and manage alignment tip labels
###############################################################################

#its alignment upload
its<-read.dna("perl_ITS_27x_20190725_mafft_ITSx.phy") #import multiple sequence alignment


# First, sample from big annotation table the rows that corespond to the sequences
# in the given alignment 
its.lab<-rownames(its)
its.annot.filter = annot.full[which(annot.full$Specimen_ID %in% its.lab), ]
its.annot.filter = drop.levels(its.annot.filter) # drop unused, "ghost" levels


# Ordering species labels by the order of sequences in alignment/distance matrix
its.annot.match <- its.annot.filter[match(its.lab, its.annot.filter$Specimen_ID),]


#Check that labels in alignment and in table annot.match are in the same order
identical(rownames(its), as.character(its.annot.match$Specimen_ID))


#Replace original tip labels in alignment by the ones needed for publication
rownames(its)<-as.character(its.annot.match$ITS_tip_labels)



###############################################################################
# Export tree & alignment with new labels for Mesquite and manuscript
###############################################################################
write.nexus.data(its, format = "dna", file= "ITS_ali_perl_Mesquite.nex", 
                 datablock=F, interleaved=F)
write.nexus(tree, file="ITS_MLtree_perl_Mesquite.nex")



###############################################################################
# Visualize the tree with annotations as boxplots
###############################################################################

# Generate treedata object
p0<-ggtree(tree)


# Adjust tree tips style:
# Replace all underscores by one blank space
p0$data$label <-gsub("_", " ", p0$data$label)


# Adjust alignment labels style:
# Replace all underscores by one blank space
rownames(its) <-gsub("_", " ", rownames(its))


# Treat the branch supports first outside of tree object
# I will plot the values higjher than 0.63 to eliminate the 
# small values (around 0) at the tip nodes. Then the figure is cleaner.
d <- p0$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
# NB! There will be a warning message "NAs introduced by coercion" because
# there were two non-numeric characters in the original tree data, namely
# the word "Root" and one empty value. 
# But this does not make a problem for our visualization.
d$label <- round(d$label, 2)
d <- d[d$label > 0.63,]


# Plot the tree and the alignment next to it
p1 <- p0 + 
  geom_tree(size=0.9) +
  geom_tiplab(size=2, align=F, hjust=-0.01)+
  geom_label_repel(data=d, aes(label=label), size = 3, nudge_x=-0.01, nudge_y=0.19, 
                   label.r = 0.05, label.padding=0.15)+
  geom_treescale(x=0.01, y=25, fontsize=3, linesize=1, width=0.05)

p1
# When plotting p1 itself and via msaplot(), 
# there will be a warning message 
# "Removed 2 rows containing missing values (geom_label_repel).
# These is due to 2 NA values introduced above. 
# It is not a problem for our visualization.

msaplot(p1, its, offset=0.52, width=3, height=0.5)
