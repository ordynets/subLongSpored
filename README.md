# subLongSpored
R code and input data files to reproduce the results of the long-spored _Subulicystidium_ paper

Each folder folders contains the R code to execute a particular task, and the input data file. 
Occasionally, output data file and graphical result are also provided. 

This code allows to assign custom, including very complex, names to the DNA sequences in the alignment and to the tips of phylogenetic tree.
This code also enables automated coloring of the tips of the phylogenetic tree according to some attribute (e.g. species name). The palette I used was successful to select 21 colors that were well distinguishable on the plotted tree!

Furthermore, this code shows how to plot the branch support values from two phylogenetic analyses onto single phylogenetic tree. The code is based on the R code from Samuel Crane (https://github.com/samuelcrane/label-node-support/) and was adjusted for plotting of the own tree with ca. 100 DNA sequences. 

This code is an easier reusable version of the code in protocol Ordynets et al. 2018 https://www.protocols.io/view/plotting-phylogenetic-tree-with-branch-supports-fr-n9fdh3n. 
