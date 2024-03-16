# This script is to make overview heatmap
# 23.11.18 Gilbert
#== environment-----------------
setwd("../../limb/3.9_wt_integrate/")
library(ComplexHeatmap)
library(RColorBrewer)

#== functions----------------------------
makeAnno <- function(metaList,colPal){
  #metaList: A list of color, list
  #colPal: color palate form RColorBrewer, str
  # example:color_organ <- makeAnno(metaSample$Organ,"Set1")
  #         color_tissue <- makeAnno(metaSample$Organ,"Set2")
  mycol=colorRampPalette(brewer.pal(7,colPal))(length(levels(metaList)))
  names(mycol) <- levels(metaList)
  return(mycol)
}
#== data prepare------------------------------
# meta
meta <- full_seurat@meta.data
showColumns <- c("Sample", "Project", "Organ", "Tissue", "Tissue.Specific.", 
             "Stage", "Gene.type", "Treatment", "Age", "Age.In.Detail.", "Machine", 
             "Species",  "Bone.Forming.Methods")
meta <- meta[,showColumns]
metaSample <- meta %>%
  group_by(Sample) %>%
  summarise(across(everything(), first, .names = "{.col}"))

# count
metaSample <- metaSample%>%column_to_rownames("Sample")
metaCount <- table(meta$Sample)%>%as.data.frame()
colnames(metaCount) <- c("Sample","Cell Count")
metaCount <- metaCount%>% column_to_rownames("Sample")
#metaCount <- log(metaCount)
levels(metaSample$Age) <- c('Organogenesis stage','Fetal stage', 'Postnatal', 'Young Adult', 'Adult','Old')


#== plot-------------------
color_project <- makeAnno(metaSample$Project,"Set3")
color_organ <- makeAnno(metaSample$Organ,"Set1")
color_tissue <- makeAnno(metaSample$Tissue,"Set2")
color_tissue_detail <- makeAnno(metaSample$Tissue.Specific.,"Set3")
color_Stage <- makeAnno(metaSample$Stage,"Pastel2")
color_Gene <- makeAnno(metaSample$Gene.type,"Pastel1")
color_treatment <- makeAnno(metaSample$Treatment,"Dark2")
color_age <- makeAnno(metaSample$Age,"Spectral")
color_ageDetail <- makeAnno(metaSample$Age.In.Detail.,"RdYlGn")
color_machine <- makeAnno(metaSample$Machine,"Accent")
color_BoneForm <- makeAnno(metaSample$Bone.Forming.Methods,"Set2")
# Annotation
ha = rowAnnotation(
  Project=metaSample$Project,
  Organ = metaSample$Organ, Age=metaSample$Age,
  AgeDetail=metaSample$Age.In.Detail.,
  Tissue=metaSample$Tissue,
  TissueDetail=metaSample$Tissue.Specific.,
  Stage=metaSample$Stage,
  Gene=metaSample$Gene.type,
  treatment=metaSample$Treatment,
  machine=metaSample$Machine,
  boneForm=metaSample$Bone.Forming.Methods,
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
    boneForm=color_BoneForm
  )
)

hm <- Heatmap(metaCount,cluster_rows = F,show_row_names = F,right_annotation = ha)
hm
haOrgan = HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE),
                       col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")))


saveRDS(hm,"processed_data/11.18_supp_overview_hm/wt_hm.Rds")
write.csv(metaSample,"processed_data/11.18_supp_overview_hm/metaSample.csv")
