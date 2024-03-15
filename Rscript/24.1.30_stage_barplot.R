unique(extendAtlas$Stage)
#extendAtlas$Stage <- extendAtlas$Stage_universal
extendAtlas$Stage <- as.character(extendAtlas$Stage)
extendAtlas$Stage[extendAtlas$Stage%in%c("Injury(5-FU)", "Injury(Non-Regeneration)", 
                                         "Injury(Radiation)", "Injury(Regeneration)")]="Injury"

#extendAtlas$Stage[is.na(extendAtlas$Stage)]="Injury"
injuryObs <- extendAtlas@meta.data[extendAtlas$C7_named=="MSC"&extendAtlas$Stage=="Injury",]
table(injuryObs$Organ)
table(injuryObs$Tissue)
injuryObs2 <- extendAtlas@meta.data[extendAtlas$update_level2=="injury MSC"&extendAtlas$Stage=="Injury",]
#injuryObs2 <- table(injuryObs2$Organ)%>%as.data.frame()
injuryObs2$Tissue <- as.character(injuryObs2$Tissue)
injuryObs2$Tissue[injuryObs2$Tissue%in%c("Femur", "Femur;Tibia", "Femurs", 
                                         "Femurs;Tibias;Ileums","Tibia")] <-  "Long bone"
tissueObs2 <- table(injuryObs2$Tissue)%>%as.data.frame()
tissueObs2=tissueObs2[tissueObs2$Freq!=0,]

#tissueObs2 <- as.character(tissueObs2$Var1)
#tissueObs2$Var1[tissueObs2$Var1%in%c("Femur", "Femur;Tibia", "Femurs", 
                                #"Femurs;Tibias;Ileums","Tibia")] <- "Long bone"
tissueObs2$Freq <- log(tissueObs2$Freq,base = 2)
tissueObs2 <- tissueObs2 %>%
  mutate(Var1 = fct_reorder(Var1, desc(Freq))) 
ggplot(tissueObs2, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +  # You can choose a different theme based on your preference
  labs(title = "Distribution of Gene Categories", y = "log2 Cell Count", x = "Tissue origin",fill="Tissue origin") +scale_fill_npg()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15,face = "bold"),  # Adjust text size here
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.title = element_text(size = 14)  # Adjust legend title size
  )
ggsave("../8.25_full_integration/result/1.1_fig6/1.30_tissueFreq.pdf",width = 6,height = 6)
ggsave("../8.25_full_integration/result/1.1_fig6/2.27_log2tissueFreq.pdf",width = 6,height = 6)
