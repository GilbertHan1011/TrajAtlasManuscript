library("TreeAndLeaf")
library("RedeR")
library("igraph")
library("RColorBrewer")
library(MicrobiotaProcess)
library(ggplot2)
library(TDbook)

library("treeio")
library("ggtree")

nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()

normaltree+ geom_nodepoint(color="red", alpha=1, size=2) 


df <- data.frame(node=c(47, 48), type=c("A", "B"))
circular_tree+ geom_nodepoint(color="red", alpha=1, size=2)+geom_tippoint(color="#FDAC4F", shape=8, size=1)+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect")

circular_tree2 =ggtree(circular_tree$data,layout = 'circular', branch.length='none',color="grey")+
  geom_nodepoint(color="#FDAC4F", alpha=1, aes(size=nodesize))+geom_tippoint(color="#FDAC4F", aes(size=nodesize))+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect",alpha=0.2)+
 
circular_tree2

meta <- wt_intergrate@meta.data
treeData <- circular_tree$data
treeDict <- data.frame(treeData$nodesize,treeData$label)
rownames(treeDict) <- treeDict$treeData.label

count1 <- table(meta$C7)%>%as.data.frame()%>%mutate(nodeSize = Freq / nrow(meta))
treeDict[count1$Var1%>%as.character(),][["treeData.nodesize"]] <- count1$nodeSize

count2 <- table(meta$C19)%>%as.data.frame()%>%mutate(nodeSize = Freq / nrow(meta))
count2 <- count2[count2$Var1%in% meta$C19,]
treeDict[count2$Var1%>%as.character(),][["treeData.nodesize"]] <- count2$nodeSize
count3 <- table(meta$C36)%>%as.data.frame()%>%mutate(nodeSize = Freq / nrow(meta))
count3 <- count3[count3$Var1%in% meta$C36,]
treeDict[count3$Var1%>%as.character(),][["treeData.nodesize"]] <- count3$nodeSize

count4 <- table(meta$C49)%>%as.data.frame()%>%mutate(nodeSize = Freq / nrow(meta))
count4 <- count4[count4$Var1%in% meta$C49,]
treeDict[count4$Var1%>%as.character(),][["treeData.nodesize"]] <- count4$nodeSize
treeDict <- is.na()

treeData$nodesize <- as.numeric(treeDict[treeData$label,]$treeData.nodesize)


circular_tree2 =ggtree(treeData,layout = 'circular', branch.length='none',color="grey")+
  geom_nodepoint(color="#FDAC4F", alpha=1, aes(size=nodesize))+geom_tippoint(color="#FDAC4F", aes(size=nodesize))+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect",alpha=0.2)+
  geom_cladelab(node=47, label="another clade",angle=-90,align=TRUE,
                geom='label', fill='lightblue',hjust='right')
circular_tree2

df <- data.frame(node=c(47, 48,49,50,51,52), type=c("MSC", "Ly6a+ MSC","Lepr+ BMSC","Fibroblast","Chondro","Pericyte"))
circular_tree2 =ggtree(treeData,layout = 'circular', branch.length='none',color="grey")+
  geom_nodepoint(color="#FDAC4F", alpha=1, aes(size=nodesize))+geom_tippoint(color="#FDAC4F", aes(size=nodesize))+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect",alpha=0.1,extend=0.5)+
  geom_cladelab(node=47, label="MSC",angle=-100,align=TRUE,
                geom='label',fill=mypal["MSC"],hjust='right',offset=0.2)+
  geom_cladelab(node=48, label="Ly6a+ MSC",angle=-100,align=TRUE,
                geom='label',fill=mypal["Ly6a+ MSC"],hjust='left',offset=0.2)+
  geom_cladelab(node=49, label="Lepr+ BMSC",angle=-100,align=TRUE,
                geom='label',fill=mypal["Lepr+ BMSC"],hjust='right',offset=0.2)+
  geom_cladelab(node=50, label="Fibroblast",angle=-100,align=TRUE,
                geom='label',fill=mypal["Fibroblast"],hjust='left',offset=0.2)+
  geom_cladelab(node=51, label="Chondro",angle=-100,align=TRUE,
                geom='label',fill=mypal["Chondro"],hjust='bottom',offset=0.2)+
  geom_cladelab(node=52, label="Pericyte",angle=-100,align=TRUE,
                geom='label',fill=mypal["Pericyte"],hjust='left',offset=0.2)+
  scale_size(range=c(0,12))+
  scale_fill_manual(values=mypal)
circular_tree2

ggsave("result/24.2.6_fig2_treeplot/tree.pdf")
circular_tree_bk <- circular_tree
