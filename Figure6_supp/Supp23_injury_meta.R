#== this script is to make supplement figure of Fig6
extendMeta <- extendAtlas@meta.data
showColumns <- c("Sample", "Project", "Core Dataset","Organ", "Tissue", "Tissue(Specific)", "Stage", "Gene type", 
                 "Treatment", "Age", "Age(In Detail)", "Machine", "Bone Forming Methods","Species")
extendMeta <- extendMeta[,showColumns]
metaSample <- extendMeta %>%
  group_by(Sample) %>%
  summarise(across(everything(), first, .names = "{.col}"))

# count
metaSample <- metaSample%>%column_to_rownames("Sample")
metaCount <- table(extendMeta$Sample)%>%as.data.frame()
metaCount=metaCount[metaCount$Freq>40,]
colnames(metaCount) <- c("Sample","Cell Count")
rownames(metaCount) <- NULL
metaCount <- metaCount%>% column_to_rownames("Sample")
metaSample <- metaSample[rownames(metaCount),]
#metaCount <- log(metaCount)
levels(metaSample$Age) <- c('Organogenesis stage','Fetal stage', 'Postnatal', 'Young Adult', 'Adult','Old')
makeAnno <- function(metaList,colPal){
  #metaList: A list of color, list
  #colPal: color palate form RColorBrewer, str
  # example:color_organ <- makeAnno(metaSample$Organ,"Set1")
  #         color_tissue <- makeAnno(metaSample$Organ,"Set2")
  mycol=colorRampPalette(brewer.pal(7,colPal))(length(unique(metaList)))
  names(mycol) <- unique(metaList)
  return(mycol)
}
#metaSample <- metaSample[rev(order(metaCount$`Cell Count`)),]
#metaCount <- metaCount[rev(order(metaCount$`Cell Count`)),]

color_project <- makeAnno(metaSample$Project,"Set3")
color_organ <- makeAnno(metaSample$Organ,"Set1")
color_tissue <- makeAnno(metaSample$Tissue,"Set2")
color_tissue_detail <- makeAnno(metaSample$`Tissue(Specific)`,"Set3")
color_Stage <- makeAnno(metaSample$Stage,"Pastel2")
color_Gene <- makeAnno(metaSample$`Gene type`,"Pastel1")
color_treatment <- makeAnno(metaSample$Treatment,"Dark2")
color_age <- makeAnno(metaSample$Age,"Spectral")
color_ageDetail <- makeAnno(metaSample$`Age(In Detail)`,"RdYlGn")
color_machine <- makeAnno(metaSample$Machine,"Accent")
color_species <- makeAnno(metaSample$Species,"Set2")
color_boneforming <- makeAnno(metaSample$`Bone Forming Methods`,"Set3")

#metaSample <- metaSample[rownames(metaCount),]
ha = rowAnnotation(
  Project=metaSample$Project,
  Organ = metaSample$Organ, Age=metaSample$Age,
  AgeDetail=metaSample$`Age(In Detail)`,
  Tissue=metaSample$Tissue,
  TissueDetail=metaSample$`Tissue(Specific)`,
  Stage=metaSample$Stage,
  Gene=metaSample$Gene.type,
  treatment=metaSample$Treatment,
  machine=metaSample$Machine,
  species=metaSample$Species,
  boneforming=metaSample$`Bone Forming Methods`,
  col = list(
    Project=color_project,
    Organ=color_organ,
    Age=color_age,
    AgeDetail=color_ageDetail,
    Tissue=color_tissue,
    TissueDetail=color_tissue_detail,
    Stage=color_Stage,
    Gene=color_Gene,
    treatment=color_treatment,
    machine=color_machine,
    species=color_species,
    boneforming=color_boneforming
  )
)


hm <- Heatmap(metaCount,cluster_rows = F,show_row_names = F,right_annotation = ha)
pdf("result/1.1_fig6/1.3_supp_meta_hm_all.pdf",width = 20,height = 12)
hm
dev.off()



meta <- read.csv()


#== heatmap overview of timp cluster----------------------------
timpMeta <- extendAtlas@meta.data[extendAtlas$update_level2=="Timp+ MSC",]
showColumns <- c("Sample", "Project", "Core Dataset","Organ", "Tissue", "Tissue(Specific)", "Stage", "Gene type", 
                 "Treatment", "Age", "Age(In Detail)", "Machine", "Species")
timpMeta <- timpMeta[,showColumns]
metaSample <- timpMeta %>%
  group_by(Sample) %>%
  summarise(across(everything(), first, .names = "{.col}"))

# count
metaSample <- metaSample%>%column_to_rownames("Sample")
metaCount <- table(timpMeta$Sample)%>%as.data.frame()
metaCount=metaCount[metaCount$Freq>40,]
colnames(metaCount) <- c("Sample","Cell Count")
rownames(metaCount) <- NULL
metaCount <- metaCount%>% column_to_rownames("Sample")
metaSample <- metaSample[rownames(metaCount),]
#metaCount <- log(metaCount)
levels(metaSample$Age) <- c('Organogenesis stage','Fetal stage', 'Postnatal', 'Young Adult', 'Adult','Old')
makeAnno <- function(metaList,colPal){
  #metaList: A list of color, list
  #colPal: color palate form RColorBrewer, str
  # example:color_organ <- makeAnno(metaSample$Organ,"Set1")
  #         color_tissue <- makeAnno(metaSample$Organ,"Set2")
  mycol=colorRampPalette(brewer.pal(7,colPal))(length(unique(metaList)))
  names(mycol) <- unique(metaList)
  return(mycol)
}
metaSample <- metaSample[rev(order(metaCount$`Cell Count`)),]
metaCount <- metaCount[rev(order(metaCount$`Cell Count`)),]

color_project <- makeAnno(metaSample$Project,"Set3")
color_organ <- makeAnno(metaSample$Organ,"Set1")
color_tissue <- makeAnno(metaSample$Tissue,"Set2")
color_tissue_detail <- makeAnno(metaSample$`Tissue(Specific)`,"Set3")
color_Stage <- makeAnno(metaSample$Stage,"Pastel2")
color_Gene <- makeAnno(metaSample$`Gene type`,"Pastel1")
color_treatment <- makeAnno(metaSample$Treatment,"Dark2")
color_age <- makeAnno(metaSample$Age,"Spectral")
color_ageDetail <- makeAnno(metaSample$`Age(In Detail)`,"RdYlGn")
color_machine <- makeAnno(metaSample$Machine,"Accent")
color_species <- makeAnno(metaSample$Species,"Set2")

#metaSample <- metaSample[rownames(metaCount),]
ha = rowAnnotation(
  Project=metaSample$Project,
  Organ = metaSample$Organ, Age=metaSample$Age,
  AgeDetail=metaSample$`Age(In Detail)`,
  Tissue=metaSample$Tissue,
  TissueDetail=metaSample$`Tissue(Specific)`,
  Stage=metaSample$Stage,
  Gene=metaSample$Gene.type,
  treatment=metaSample$Treatment,
  machine=metaSample$Machine,
  species=metaSample$Species,
  col = list(
    Project=color_project,
    Organ=color_organ,
    Age=color_age,
    AgeDetail=color_ageDetail,
    Tissue=color_tissue,
    TissueDetail=color_tissue_detail,
    Stage=color_Stage,
    Gene=color_Gene,
    treatment=color_treatment,
    machine=color_machine,
    species=color_species
  )
)


hm <- Heatmap(metaCount,cluster_rows = F,show_row_names = F,right_annotation = ha)
pdf("result/1.1_fig6/1.3_supp_meta_hm_timp1.pdf",width = 13,height = 5)
hm
dev.off()