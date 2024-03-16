unique(wt_integrate@meta.data$Age)
wtMeta <- readRDS("../important_processed_data/9.10_wt_meta.Rds")
cellCount <- table(wtMeta$Age)%>%as.data.frame()
cellCountDetail <- table(wtMeta$Age.In.Detail.)%>%as.data.frame()
cellCountDetail$Freq <- log(cellCountDetail$Freq,base = 2)
colnames(cellCountDetail) <- c("age","log2Count")

cellCountDetail$age <- factor(cellCountDetail$age,levels = c("E9", "E10", "E10.5", "E11", "E11.5", "E12", "E12.5", "E13", 
                                      "E13.5", "E14", "E14.5", "E15", "E15.5", "E16", "E16.5", "E17", 
                                      "E17.5", "E18", "E18.5", "P2", "P5", "P7.5", "P10", "P13", "P19", 
                                      "P21", "3W", "P23", "P28","4W",  "1M", "1.5M", "2M", "3M", "4M", 
                                      "Adult", "16M", "18M"))
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
ggsave("result/24.2.11_Fig2_supp/cellCount.pdf",width = 15,height = 4)    

cellCountTissue <- table(wtMeta$Tissue.Specific.)%>%as.data.frame()
cellCountTissue$Freq <- log(cellCountTissue$Freq,base = 2)
levels(cellCountTissue$Var1) <- c( "Calvarial bone", "Coronal suture","Sagittal suture",  "Frontal suture",
                                   "Lambdoid suture", "Mandible mesenchyme", "Maxillary mesenchyme",
                                  "Cranial mesenchyme",
                                  "Limb Bud",  "Perichondrial",
                                  "Bone", "Bone marrow","Diaphysis", "Diaphysis and Epiphysis", 
                                  "Epiphysis", "Growth plate", "Stroma","Metaphysis"
                                    )
colnames(cellCountTissue) <- c("tissue","log2Count")
mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(nrow(cellCountTissue))
ggplot(cellCountTissue, aes(x=tissue, y=log2Count,fill=tissue)) +
  geom_bar(stat = 'identity')+theme_classic()+
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
ggsave("result/24.2.11_Fig2_supp/cellCountTissue.pdf",width = 8,height = 4)    

projectOrgan <- unique(wtMeta[c("Sample","Organ","Age")])
c7Project <- table(wtMeta$Sample,wtMeta$C7_named)%>%as.data.frame()
colnames(c7Project) <- c("Sample","Anno","Freq")
c7Project <- left_join(c7Project,projectOrgan)

# ggplot(c7Project, aes(x=Project, y=Freq,fill=Anno)) + 
#   geom_bar(position="stack", stat="identity")+theme_classic()
# ggplot(data, aes(fill=condition, y=value, x=specie)) + 
#   geom_bar(position="stack", stat="identity")
c7Project$Organ[c7Project$Sample=="PerichondrialP21_Matsushita_Prrx1creE11.5"] <- "Limb_embryo"
c7Project$Organ <- as.character(c7Project$Organ)
c7Project$Organ[c7Project$Organ=="Limb_embryo"] <- "Limb bud"
c7Project$Organ[c7Project$Organ=="Limb_adult"] <- "Long bone"
mypal <- c("#E64B35B2", "#4DBBD5B2", "#3C5488B2","#00A087B2",  "#F39B7FB2","#e377c2ff","#00ff00ff")

names(mypal) <- c("Ob", "Fibroblast","Chondro","MSC", "Lepr+ BMSC", 
                  "Pericyte", "Ly6a+ MSC")
my_palette <- brewer.pal(n = 7, name = "Set1") 
ggplot(c7Project, aes(x=Sample, y=Freq,fill=Anno)) + 
  geom_bar(position="fill", stat="identity")+ facet_wrap(vars(Organ,Age),scales = "free")+scale_fill_manual(values=my_palette)+
  theme_classic()+labs(x = NULL)

ggplot(c7Project, aes(x = Sample, y = Freq, fill = Anno)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(Organ~ Age, scales = "free") +
  scale_fill_manual(values = mypal) +  # Set fill colors with manual scale
  theme_classic() +
  labs(x = NULL) +
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.key = element_rect(fill = "white")
  )
ggsave("result/12.19_figure2_OPC_plot/age_organ_percentage.pdf",width = 8,height =6)

ggsave("result/12.19_figure2_OPC_plot/age_organ_percentage.pdf",width = 6,height = 5)

c7Project$Age <- factor(c7Project$Age,levels = c("Organogenesis stage","Fetal stage", "Postnatal", "Young Adult", "Adult", "Old"))
ggplot(c7Project, aes(x = Sample, y = Freq, fill = Anno)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_grid(Organ~ Age, scales = "free") +
  scale_fill_manual(values = mypal) +  # Set fill colors with manual scale
  theme_classic() +
  labs(x = NULL) +
  theme(
    plot.title = element_text(size = 20),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.key = element_rect(fill = "white")
  )
ggsave("result/12.19_figure2_OPC_plot/age_organ_percentage_grid.pdf",width = 10,height = 4)

write.csv(c7Project,"processed_data/12.23_fig2_supp/figure2d_facet_grid.csv")





