
scList <- list(ls(pattern = "sc_+"))
# Create an empty list to store the objects
objList <- list()

# Iterate over the character strings and retrieve the objects
for (i in 1:length(scList[[1]])) {
  objList[[i]] <- get(scList[[1]][i])
}

# Print the list of objects
objList


anno <- readxl::read_xlsx("data/ssc_anno.xlsx")
names(anno$`Stem cell name`)

new.id <- c( "Acta2+ BMSC", "Adipoq+ BMSC", "Axin2+ MSC","Cdh2+ BMSC",  "Clec11a+ BMSC","Hypertrophic chondrocyte", 
             "Chondrocyte","Cxcl12+ BMSC", "Early perichondrial cell", "Fox2a+ chondrocyte",
             "Gli1+ MSC", "Grem1+ BMSC", "Hhip+ MSC",  "Cd168+ skeletal stem/progenitor cell", "Lepr+ BMSC","Mcam+ BMSC",  
             "Msx2+ MSC", "Pdgfra+ BMSC","Prrx1+ MSC","Pthlh+ chondrocyte", "pvSSC",  "ocSSC")
names(objList) <- new.id
saveRDS(objList,"processed_data/4.22_highlightcell.Rds")
objList <- readRDS("processed_data/4.22_highlightcell.Rds")

DimPlot(full_seurat,cells.highlight = objList["Axin2+ MSC"])
rm(limbAdult)
rm(limbEmbryo)
rm(limbEmbryo2)
test <- full_seurat
test <- CreateSeuratObject(counts = t(full_seurat@reductions$X_scANVI@cell.embeddings))
test@meta.data <- full_seurat@meta.data
test@reductions <- full_seurat@reductions
DimPlot(test)
saveRDS(test,"data/mini_seurat.Rds")


objList <- readRDS("processed_data/4.22_highlightcell.Rds")
sc_hes1_cell <- list(sc_hes1_cell)
names(sc_hes1_cell) <- "Hes1+ Perichondrium"
objList <- append(objList,sc_hes1_cell)
cellArray <- unlist(objList)
cell_freq <- table(cellArray)
cell_freq <- as.data.frame(cell_freq)
sscFreq <- table(cell_freq$Freq)
sscFreq <- as.data.frame(sscFreq)
hist(sscFreq)
names(sscFreq) <- c("Freq","Cell number")
# Create a barplot using ggplot2
ggplot(sscFreq, aes(x = Freq, y = `Cell number`, fill = Freq)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "SSC Freq", x = "Freq", y = "Cell number")+
  theme_minimal() +
  scale_fill_manual(values  =  rev(brewer.pal(n = 6, name = "RdBu")))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), # Make title bold
        axis.text =  element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"), # Make x-axis label bold
        axis.title.y = element_text(face = "bold"))
ggsave("result/4.22_stem_cell_anno/stemcell_freq.pdf")
cellList <- list()
for (i in 1:6){
  cellList[[i]] <- as.character(cell_freq$cellArray[cell_freq$Freq==i])
  DimPlot(full_seurat,cells.highlight = cellList[[i]])+labs(title = paste0("Stem cells that have been studied " ,i, " times"))
  ggsave(paste0("result/4.22_stem_cell_anno/anno_",i, "_time_ssc.pdf"))
}

library(VennDiagram)

venn.diagram(
  objList,
  #category.names = c("Fruits", "Fruits and Grapes", "Fruits and Kiwis"),
  filename = NULL,
  # output = FALSE,
  # fill = c("cornflowerblue", "green", "yellow"),
  # alpha = rep(0.5, 3),
  # label.col = rep("white", 3),
  # cex = 2.5,
  # fontfamily = "serif"
)
library(UpSetR)
UpSetR::upset(fromList(objList),nsets = 23,order.by = "freq", sets.bar.color = "#69b3a2",text.scale=2, point.size = 3.5, line.size = 2,main.bar.color = "#69b3a2")
saveRDS(objList,"processed_data/4.22_highlightcell.Rds")
