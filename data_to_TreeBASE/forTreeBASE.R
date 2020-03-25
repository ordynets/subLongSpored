###############################################################################
# Load required packages
###############################################################################
library(here)
library(ape)
library(Hmisc)
library(gdata)
library(ape)
library(geiger)
#library(phyloch) # Can be ignored if using the conformat=simple option in MrBayes sumt command

# Set the current directory as working directory
set_here()

# Read in metadata 
annot.full <- read.csv ("SUB_data_processing_20190731.csv",  sep="\t")



###############################################################################
# Read in Bayes tree and adjust the tip labels
###############################################################################

# Read in Bayes tree
read.nexus("SUB_BY_outp.con.tre")->bayesTree2 #Reads in the .con file that results from analyses in MrBayes.
bayesTree2[[1]]->bayesTree.unr #Extracts one of the two trees in the .con file.
bayesTree.unr$node.label<-round(as.numeric(bayesTree.unr$node.label), digits=2) # round pp values

#This is used to see node labels and used the particular node to reroot around it (see below)
plot(bayesTree.unr, show.node.label = T) 
nodelabels()


# Proper rooting of the Bayes tree
bayesTree<-root(bayesTree.unr, "KAS_L_1860", resolve.root = T) #NB! Check node in each analysis! 
#plot(bayesTree)
is.rooted(bayesTree)

# First, sample from big annotation table the rows that correspond to the tips in the given tree
by.tip.lab<-bayesTree$tip.label
annot.filter.by <- annot.full[which(annot.full$Specimen_ID %in% by.tip.lab), ]
annot.filter.by <- drop.levels(annot.filter.by) # drop unused, "ghost" levels


# Ordering species labels by the order of sequences in alignment/distance matrix
annot.match.by <- annot.filter.by[match(by.tip.lab, annot.filter.by$Specimen_ID),]


#Check that labels in tree and in table annot.match are in the same order
identical(bayesTree$tip.label, as.character(annot.match.by$Specimen_ID))


#Replace original tip labels in Bayes tree by the ones needed for publication
bayesTree$tip.label<-as.character(annot.match.by$CONCAT_tip_labels)


###############################################################################
# Read in RAxML tree and adjust the tip labels
###############################################################################

# Read in RAxML tree
read.tree("PhyML_newick_tree.nhx")->bootTree
bootTree$node.label<-round(as.numeric(bootTree$node.label), digits=2) # round aLRT values


# First, sample from big annotation table the rows that corespond to the tips in the given tree
ml.tip.lab<-bootTree$tip.label
annot.filter.ml <- annot.full[which(annot.full$Specimen_ID %in% ml.tip.lab), ]
annot.filter.ml <- drop.levels(annot.filter.ml) # drop unused, "ghost" levels


# Ordering species labels by the order of sequences in alignment/distance matrix
annot.match.ml <- annot.filter.ml[match(ml.tip.lab, annot.filter.ml$Specimen_ID),]

#Check that labels in tree and in table annot.match are in the same order
identical(bootTree$tip.label, as.character(annot.match.ml$Specimen_ID))

#Replace original tip labels in Bayes tree by the ones needed for publication
bootTree$tip.label<-as.character(annot.match.ml$CONCAT_tip_labels)


###############################################################################
# # Load alignment and adjust sequence labels for TreeBASE
###############################################################################

#concat alignment upload
concat<-read.dna("CONCAT_all_20190725_v3_noQuestionmarks.phy") #import multiple sequence alignment


# First, sample from big annotation table the rows that corespond to the sequences
# in the given alignment 
concat.lab<-rownames(concat)
concat.annot.filter = annot.full[which(annot.full$Specimen_ID %in% concat.lab), ]
concat.annot.filter = drop.levels(concat.annot.filter) # drop unused, "ghost" levels


# Ordering species labels by the order of sequences in alignment/distance matrix
concat.annot.match <- concat.annot.filter[match(concat.lab, concat.annot.filter$Specimen_ID),]


#Check that labels in alignment and in table annot.match are in the same order
identical(rownames(concat), as.character(concat.annot.match$Specimen_ID))


#Replace original tip labels in alignment by the ones needed for publication
rownames(concat)<-as.character(concat.annot.match$CONCAT_tip_labels)

# Adjust alignmen labels style for TreeBASE:
# Replace all underscores and colons by underscore
rownames(concat) <-gsub(" ", "_", rownames(concat))
rownames(concat) <-gsub(":", "_", rownames(concat))
rownames(concat) <-gsub("-", "_", rownames(concat))


###############################################################################
# Adjust Tree tip labels for TreeBASE
###############################################################################
bayesTree$tip.label <-gsub(" ", "_", bayesTree$tip.label)
bayesTree$tip.label <-gsub(":", "_", bayesTree$tip.label)
bayesTree$tip.label <-gsub("-", "_", bayesTree$tip.label)

bootTree$tip.label <-gsub(" ", "_", bootTree$tip.label)
bootTree$tip.label <-gsub(":", "_", bootTree$tip.label)
bootTree$tip.label <-gsub("-", "_", bootTree$tip.label)

###############################################################################
# Export tree & alignment with labels for TreeBASE
###############################################################################
write.nexus.data(concat, format = "dna", file= "CONCAT_ali_all_Mesquite.nex", 
                 datablock=F, interleaved=F)
write.nexus(bayesTree, file="CONCAT_BYtree_all_Mesquite.nex")
write.nexus(bootTree, file="CONCAT_MLtree_all_Mesquite.nex")

