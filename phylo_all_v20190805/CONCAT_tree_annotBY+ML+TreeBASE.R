###############################################################################
# Load required packages
###############################################################################
library(here)
library(ape)
library(Hmisc)
library(gdata)
library(pals)
library(ape)
library(geiger)
#library(phyloch) # Can be ignored if using the conformat=simple option in MrBayes sumt command

# Set the current directory as working directory
set_here()

# Read in data 
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
# Generate color palette for annotation of tips on Bayesian tree topology
###############################################################################

# Create palette of distinct colors
# the palette does not change between runs
mycol <- alphabet2(n = length(levels(as.factor(annot.match.by$Species_name_v2))))



###############################################################################
# Defining functions which will confront clades of one tree to 
# clades of another tree
###############################################################################
# copied from 
# https://github.com/samuelcrane/label-node-support/blob/cf4671dc01eab90e9c7d06aae29fabae2df0f834/labelNodeSupport.r
# labelNodeSupport: map node support values from parsimony, likelihood, 
# and Bayesian phylogenies onto a single target tree.
# Samuel Crane wrote:
# "This script was put together by Rich Glor and first published (as far as I can tell) in 2008 
# on the Dechronization blog # http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html
# I've updated the code (some of the syntax was depreciated), made some corrections, 
# and am in the process of expanding the function by
# (a) adding in the ability to map a third set of support values and
# (b) deal with the verbose (as well as the simple) output of MrBayes sumt consensus command [SNC]"


# The getAllSubTrees function below is a necessary subfunction that atomizes a tree into each 
# individual subclade and was provided compliments of Luke Harmon. [RG]
getAllSubtrees <- function(phy, minSize=2) 
{
  res <- list() 
  count = 1 
  ntip <- length(phy$tip.label) 
  for(i in 1:phy$Nnode) 
  { 
    l <- tips(phy, ntip+i) 
    bt <- match(phy$tip.label, l) 
    if(sum(is.na(bt)) == 0) 
    {
      st <- phy 
    } 
    else st <- drop.tip(phy, phy$tip.label[is.na(bt)]) 
    if(length(st$tip.label)>=minSize) 
    { 
      res[[count]] <- st 
      count <- count+1 
    }
  } 
  res
}


# The plotBayesBoot function below plots both posterior probability and bootstrap values on each 
# node of the consensus tree obtained from your Bayesian analysis. Bootstrap values will appear in 
# bold text immediately below and to the left of the node they support, whereas Bayesian posterior 
# probabilies will appear in regular face above and to the left of the node. [RG]

plotBayesBoot <- function(bayesTree,bootTree) 
{
  getAllSubtrees(bayesTree) -> bayesSub
  getAllSubtrees(bootTree) -> bootSub
  bootList <- matrix("NA", Nnode(bayesTree), 1)
  
  #The commands below compare all the subclades in the Bayes tree to all the subclades 
  #in the bootstrap tree, and vice versa, and identifies all those clades that are identical.
  for(i in 1:Nnode(bayesTree)) 
  {
    for(j in 1:Nnode(bootTree)) 
    {
      match(bayesSub[[i]]$tip.label[order(bayesSub[[i]]$tip.label)], bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)]) -> shared
      match(bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)], bayesSub[[i]]$tip.label[order(bayesSub[[i]]$tip.label)]) -> shared2
      if(sum(is.na(c(shared,shared2)))==0) 
      {
        bootTree$node.label[j] -> bootList[i]
      }
    }
  }
  plot(ladderize(bayesTree, right=F), cex=0.4, font=1, edge.width = 1.7, use.edge.length=TRUE, direction='r', 
       tip.color=mycol[as.factor(annot.match.by$Species_name_v2)]) #Plots your Bayesian consensus tree
  nodelabels(bayesTree$node.label, adj=c(1.3, -0.6), frame="n", cex=0.5, font=1) #Adds posterior probability values to the tree. Change the 'cex' value to make the fond smaller or larger. A value of 1 will give you a readable result in the R quartz window, but a value closer to 0.25 might be better for publication)
  nodelabels(bootList, adj=c(1.3, 1.4), frame="n", cex=0.5, font=1) #Adds bootstrap values
  
  legend(0.000, 100, legend=levels(as.factor(annot.match.by$Species_name_v2)), col=mycol, 
         pch=15, pt.cex=1.2, bty="n", cex=0.4,  title=expression(bold("Subulicystidium species")), title.adj = 1.4, 
         y.intersp=0.9)
  #par(op) #reset graphical parameters to defaults
  
}



###############################################################################
# Final picture
###############################################################################


par(mar=c(0,0,0,0))
plotBayesBoot(bayesTree, bootTree) # For two trees
add.scale.bar(0.01, 55, length=0.01, lwd=2, cex=0.7)
dev.off()

# End