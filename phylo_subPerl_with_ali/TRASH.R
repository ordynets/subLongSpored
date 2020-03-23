
#creat a tree as an R object 
#p1<-ggtree(tree, aes(x, y)) + 
#geom_tiplab(size=1, align=F)+
#geom_text(aes(x=branch,label=label), vjust=-0.5) # this to activate when the
#file with branch support values will be available
#geom_treescale(x=0.01, y=27, fontsize=3)


tree$node.label<-as.numeric(tree$node.label)
tree$node.label<-round(tree$node.label, 2)


p1<-ggtree(tree, size=1) + 
  geom_tiplab(size=3, align=F, hjust=-0.01)+
  geom_text2(aes(subset = !isTip, label=label), size=2, vjust=-0.5, hjust=1.2)


msaplot(p1, its, offset=0.5, width=3, height=0.6)

