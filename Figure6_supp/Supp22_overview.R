extendMeta <- readRDS("../important_processed_data/9.10_merge_meta.Rds")
unique(extendMeta$Age.In.Detail.)
extendAtlas <- readRDS("../important_processed_data/12.30_full_latent.Rds")
extendAge <- unique(extendAtlas$`Age(In Detail)`)
wtAge <- unique(wt_integrate$Age.In.Detail.)
extendAtlas$`Age(In Detail)` <- factor(extendAtlas$`Age(In Detail)`,
                                       levels = c("5 WPC","8 WPC","E9", "E10.5", "E10","12 WPC", "13 WPC", "14 WPC", "15 WPC", "16 WPC",
                                                  "17 WPC", "19 WPC", "E11", "E11.5", "E12", "E12.5", "E13", 
                                                  "E13.5", "E14", "E14.5", "E15", "E15.5", "E16", "E16.5", "E17", 
                                                  "E17.5", "E18", "E18.5", "P2","P3",  "P4-P5","P5","P6","P7", "P7.5",
                                                  "P10","P11", "P13","P14","P18", "3W", 
                                                  "P19", "P21", "P23","P25", "P28", "2-4M","4W", "1M","5W", "6W","P42","1.5M","7W","8W",
                                                  "P56","2M","8-12W", "9-12W", "10-12W","11-12W", "10-14W",
                                                  "10W","11W","12W","3M","14W", "4M", 
                                                  "Adult", "12M","12-14M","16M", "18M"))
cellCountDetail <- table(extendAtlas$`Age(In Detail)`)%>%as.data.frame()
cellCountDetail$Freq <- log(cellCountDetail$Freq,base = 2)
colnames(cellCountDetail) <- c("age","log2Count")
mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(nrow(cellCountDetail))
ggplot(cellCountDetail, aes(x=age, y=log2Count,fill=age)) +
  geom_bar(stat = 'identity')+theme_classic()+scale_fill_manual(values=rev(mycolor))+
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 15, color = "black",angle = 90),
    axis.text.y = element_text(size = 15, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 15),
    legend.key = element_rect(fill = "white"))
ggsave("result/24.2.26_fig6_supp//cellCount_age.pdf",width = 18,height = 4)



tissueMeta <- unique(extendAtlas@meta.data[c("Organ","Tissue(Specific)")])
tissueMeta <- tissueMeta%>%
  arrange(Organ)
extendAtlas$`Tissue(Specific)` <- factor(extendAtlas$`Tissue(Specific)`,
                                         levels =c("Digit bone", "Suture mesenchyme", "Cranial mesenchyme", "Calvarial bone", 
                                                   "Mandible mesenchyme", "Maxillary mesenchyme", "Lambdoid suture", 
                                                   "Periodontium of Mandibular molars", "Parietal bone", 
                                                   "Diaphysis", "Diaphysis and Epiphysis", "Bone", "Stroma", "Epiphysis", 
                                                   "Metaphysis", "Growth plate",
                                                   "Metaphysis and Diaphysis", "Femurs;Tibias;Ileums", "Bone marrow", 
                                                   "Tendon enthesis", "Periosteal",  "Knee Synovial Joint", 
                                                   "Synovial Joint", "Skeletal Muscle", "Limb Bud", 
                                                   "Perichondrial","Repair Callus of Rib Bones", 
                                                   "Subdermal ectopic mass cells"))

cellCountDetail <- table(extendAtlas$`Tissue(Specific)`)%>%as.data.frame()
cellCountDetail$Freq <- log(cellCountDetail$Freq,base = 2)
colnames(cellCountDetail) <- c("tissue","log2Count")

ggplot(cellCountDetail, aes(x=tissue, y=log2Count,fill=tissue)) +
  geom_bar(stat = 'identity')+theme_classic()+
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10, color = "black",angle = 90),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = "white"))
ggsave("result/24.2.26_fig6_supp//cellCount_tissue.pdf",width = 10,height =6)



tissueMeta <- unique(extendAtlas@meta.data[c("Organ","Tissue(Specific)")])
tissueMeta <- tissueMeta%>%
  arrange(Organ)
extendAtlas$`Tissue(Specific)` <- factor(extendAtlas$`Tissue(Specific)`,
                                         levels =c("Digit bone", "Suture mesenchyme", "Cranial mesenchyme", "Calvarial bone", 
                                                   "Mandible mesenchyme", "Maxillary mesenchyme", "Lambdoid suture", 
                                                   "Periodontium of Mandibular molars", "Parietal bone", 
                                                   "Diaphysis", "Diaphysis and Epiphysis", "Bone", "Stroma", "Epiphysis", 
                                                   "Metaphysis", "Growth plate",
                                                   "Metaphysis and Diaphysis", "Femurs;Tibias;Ileums", "Bone marrow", 
                                                   "Tendon enthesis", "Periosteal",  "Knee Synovial Joint", 
                                                   "Synovial Joint", "Skeletal Muscle", "Limb Bud", 
                                                   "Perichondrial","Repair Callus of Rib Bones", 
                                                   "Subdermal ectopic mass cells"))

cellCountDetail <- table(extendAtlas$Species)%>%as.data.frame()
cellCountDetail$Freq <- log(cellCountDetail$Freq,base = 2)
colnames(cellCountDetail) <- c("Species","log2Count")
mycolor<-colorRampPalette(brewer.pal(8,'Set1'))(nrow(cellCountDetail))
ggplot(cellCountDetail, aes(x=Species, y=log2Count,fill=Species)) +
  geom_bar(stat = 'identity')+theme_classic()+scale_fill_manual(values=rev(mycolor))+
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black",angle = 90),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = "white"))
ggsave("result/24.2.26_fig6_supp/cellCount_species.pdf",width = 5,height =4)

cellCountDetail <- table(extendAtlas$Stage)%>%as.data.frame()
cellCountDetail$Freq <- log(cellCountDetail$Freq,base = 2)
colnames(cellCountDetail) <- c("Stage","log2Count")
mycolor<-colorRampPalette(brewer.pal(8,'Set1'))(nrow(cellCountDetail))
ggplot(cellCountDetail, aes(x=Stage, y=log2Count,fill=Stage)) +
  geom_bar(stat = 'identity')+theme_classic()+scale_fill_manual(values=rev(mycolor))+
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black",angle = 90),
    axis.text.y = element_text(size = 20, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = "white"))
ggsave("result/24.2.26_fig6_supp/cellCount_treat.pdf",width = 10,height =5)


cellCountDetail <- table(extendAtlas$`Gene type`)%>%as.data.frame()
cellCountDetail$Freq <- log(cellCountDetail$Freq,base = 2)
colnames(cellCountDetail) <- c("geneType","log2Count")
mycolor<-colorRampPalette(brewer.pal(8,'Set1'))(nrow(cellCountDetail))
ggplot(cellCountDetail, aes(x=geneType, y=log2Count,fill=geneType)) +
  geom_bar(stat = 'identity')+theme_classic()+scale_fill_manual(values=rev(mycolor))+
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black",angle = 90),
    axis.text.y = element_text(size = 20, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 20),
    legend.key = element_rect(fill = "white"))
ggsave("result/24.2.26_fig6_supp/cellCount_geneType.pdf",width = 8,height =5)


cellCountDetail <- table(extendAtlas$Machine)%>%as.data.frame()
cellCountDetail$Freq <- log(cellCountDetail$Freq,base = 2)
colnames(cellCountDetail) <- c("Machine","log2Count")
mycolor<-colorRampPalette(brewer.pal(8,'Set1'))(nrow(cellCountDetail))
ggplot(cellCountDetail, aes(x=Machine, y=log2Count,fill=Machine)) +
  geom_bar(stat = 'identity')+theme_classic()+scale_fill_manual(values=rev(mycolor))+
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black",angle = 90),
    axis.text.y = element_text(size = 20, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = "white"))
ggsave("result/24.2.26_fig6_supp/cellCount_Machine.pdf",width = 8,height =5)


#== featureplot------------


DimPlot(extendAtlas,group.by = "Stage")


