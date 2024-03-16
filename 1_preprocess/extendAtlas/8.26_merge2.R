metadata <- readxl::read_xlsx("data/8.18_sampleMeta.xlsx")
setdiff(metadata$Sample,mergeAll$orig.ident)
unique(mergeAll$orig.ident)

mergeAll$orig.ident[mergeAll$orig.ident=="Frontal_Holmes_Fgfr2E18.6"] <- "Frontal_Holmes_Fgfr2E18.5"
intersect(colnames(Bmsc2019_Regev_Leu),colnames(mergeAll))
length(mergeAll$Condition==c("1:ctrl",  "4:ctrl","9:ctrl","5:ctrl", 
                             "15:ctrl",  "8:ctrl", "3:ctrl", 
                             "13:ctrl", "16:ctrl","7:ctrl", 
                             "2:ctrl", "12:ctrl","10:ctrl", "17:ctrl"))
length(mergeAll$Condition[mergeAll$Condition==c("1:ctrl",  "4:ctrl","9:ctrl","5:ctrl", 
                             "15:ctrl",  "8:ctrl", "3:ctrl", 
                             "13:ctrl", "16:ctrl","7:ctrl", 
                             "2:ctrl", "12:ctrl","10:ctrl", "17:ctrl")])
leuName <- na.omit(colnames(mergeAll)[mergeAll$Condition=="ctrl"])

mergeAll$orig.ident[leuName] <- gsub("samp", "LeuCon", mergeAll$orig.ident[leuName])

Septoclasts_Kishor <- readRDS("../data/cleandata/Septoclasts_Kishor_Fracture.Rds")
Septoclasts_Kishor <- DietSeurat(Septoclasts_Kishor)
unique(Septoclasts_Kishor$orig.ident)
new.id <- c("Septoclasts_Kishor_Con1",
            "Septoclasts_Kishor_Con2",
            "Septoclasts_Kishor_Con3",
            "Septoclasts_Kishor_Con4",
            "Septoclasts_Kishor_Fracture1",
            "Septoclasts_Kishor_Fracture2",
            "Septoclasts_Kishor_Fracture3",
            "Septoclasts_Kishor_Fracture4"
)
names(new.id) <- unique(Septoclasts_Kishor$orig.ident)
Idents(Septoclasts_Kishor) <- Septoclasts_Kishor$orig.ident
Septoclasts_Kishor <- RenameIdents(Septoclasts_Kishor,new.id)
Septoclasts_Kishor$orig.ident <- Idents(Septoclasts_Kishor)


LimbG610C_Gorrell <- readRDS("../data/cleandata/LimbG610C_Gorrell_calvaria.Rds")
unique(LimbG610C_Gorrell$orig.ident)
LimbG610C_Gorrell <- DietSeurat(LimbG610C_Gorrell)


digitColnames <- colnames(mergeAll[,mergeAll$orig.ident%in%c("1","2","3","4","5","6","7","8","9","10")])
digitColnames<-  sub("_[^_]*$", "", digitOrigName)
digitColnames<-  sub("_[^_]*$", "", digitColnames)
digitOrigName <-  colnames(mergeAll[,mergeAll$orig.ident%in%c("1","2","3","4","5","6","7","8","9","10")])
digitGenName <- digitOrigName[digitColnames%in%colnames(Digit2019_Storer_Gen)]
digitNonGenName <- digitOrigName[digitColnames%in%colnames(Digit2019_Storer_NonGen)]

digitColnames%in%colnames(Digit2019_Storer_Gen)
digitCols <- colnames(mergeAll[,mergeAll$orig.ident%in%c("1","2","3","4","5","6","7","8","9","10")])
names(digitCols) <- digitColnames


mergeAll$orig.ident[digitGenName] <- Digit2019_Storer_Gen$orig.ident
mergeAll$orig.ident[digitNonGenName] <- as.character(Digit2019_Storer_NonGen$orig.ident)
mergeAll <- merge(mergeAll,y = list(Septoclasts_Kishor,LimbG610C_Gorrell))

source("../function/seurat_utils.R")


#== create metadata-------------
mergeAll$Sample <- mergeAll$orig.ident


for (i in colnames(metadata[2:ncol(metadata)])){
  mergeAll <- createSampleMeta(mergeAll,metadata,i)
}
createSampleMeta$short_id <- createSampleMeta$Project%>%
  stringr::str_split("_")%>%
  map(.,function(x) paste(substr(x,1,2),collapse = ""))%>%unlist()
# osteoMerge <- createSampleMeta(osteoMerge,osteoMeta,"short_id")
# 
unique(mergeAll$Project)
dir.create("../unimportant_processed_data/8.25_integrateAll")
saveRDS(mergeAll@meta.data,"../unimportant_processed_data/8.25_integrateAll/8.25_metaBk.Rds")
selectMeta <- c(colnames(metadata),"orig.ident", "nCount_RNA", "nFeature_RNA",
                "paper_label", "coarse_label",  "scDblFinder_class","short_id")
# tset <- osteoMeta$Project%>%
#   stringr::str_split("_")%>%
#   map(.,function(x) paste(substr(x,1,1),collapse = ""))

# osteoMerge$short_id <- osteoMerge$Project%>%
#   stringr::str_split("_")%>%
#   map(.,function(x) paste(substr(x,1,1),collapse = ""))
mergeAll@meta.data <- mergeAll@meta.data[,selectMeta]
saveRDS(mergeAll@meta.data,"../unimportant_processed_data/8.25_integrateAll/8.25_metaBk.Rds")
saveRDS(mergeAll,"../important_processed_data/8.27_merge_78w.Rds")
library(SeuratDisk)




test1 <- mergeAll[,mergeAll$Sample=="Suture2021_Farmer_E17"]
test2 <- mergeAll[,mergeAll$Sample=="SutureP0_Tower_WT"]
test3 <- mergeAll[,mergeAll$Sample=="LimbSSC_Chan"]


test1@assays$RNA@counts <- test1@assays$originalexp@counts


origin1 <- c("Suture2021_Farmer_E17", "Suture2021_Farmer_E15", "CranioSoxc_Angelozzi_WTE11.5", 
            "CranioSoxc_Angelozzi_WTE12.5", "CranioSoxc_Angelozzi_WTE13.5", 
            "CranioSoxc_Angelozzi_WTE15.5", "CranioSoxc_Angelozzi_WTE17.5", 
            "CalvariaP4_Ayturk", "Mandible2020_Chai", "Frontal_Holmes_PN10_Anterior", 
            "Frontal_Holmes_WTE16.5_1", "Frontal_Holmes_WTE16.5_2", "Frontal_Holmes_WTE18.5_1", 
            "Frontal_Holmes_WTP28", "Frontal_Holmes_PN10_Posterior", "Frontal_Holmes_WTE18.5_2", 
            "Frontal_Holmes_WTE18.5_3", "Frontal_Holmes_WTE18.5_4", "Maxillary_Bian_E10.5", 
            "Maxillary_Bian_E11.5", "Maxillary_Bian_E12.5", "Maxillary_Bian_E14.5", 
            "Mesenchymal2022_Zhang_E9", "Mesenchymal2022_Zhang_E10", "Mesenchymal2022_Zhang_E11.5", 
            "Mesenchymal2022_Zhang_E14.5", "lambdoid_Holmes_E18", "lambdoid_Holmes_P10", 
            "lambdoid_Holmes_P28", "sagittal_Holmes_E18", "sagittal_Holmes_P10", 
            "sagittal_Holmes_P28", "coronal_Holmes_E16_1", "coronal_Holmes_E16_2", 
            "coronal_Holmes_E18_1", "coronal_Holmes_E18_2", "coronal_Holmes_P10", 
            "coronal_Holmes_P28", "Mandible2020_Chai_E10", "Mandible2020_Chai_E12", 
            "ChondroOsteo_Long", "Forelimb_He_E13.5", "Forelimb_He_E11", 
            "Forelimb_He_E13_1", "Forelimb_He_E10.5", "Forelimb_He_E13_2", 
            "Forelimb_He_E15_1", "Forelimb_He_E14", "Forelimb_He_E15_2", 
            "Forelimb_He_E12", "LimbG610C_Gorrell_femurWT2", "LimbG610C_Gorrell_femurWT3", 
            "LimbG610C_Gorrell_femurWT4", "LimbG610C_Gorrell_femurWT1", "LimbMouse2019_Kelly_E11", 
            "LimbMouse2019_Kelly_E13", "LimbMouse2019_Kelly_E15", "LimbMouse2019_Kelly_E18", 
            "PerichondrialE13.5_Matsuacshita", "PerichondrialP21_Matsushita_Prrx1creE11.5", 
            "Ablation_Matsushita_abl14", "Ablation_Matsushita_abl7con1", 
            "Ablation_Matsushita_abl7con2", "Ablation_Matsushita_treat1", 
            "Ablation_Matsushita_treat2", "Ablation_Matsushita_cxcl1", "Ablation_Matsushita_cxcl2", 
            "Bmsc2019_Regev_samp1", "Bmsc2019_Regev_samp2", "Bmsc2019_Regev_samp3", 
            "Bmsc2019_Regev_samp4", "Bmsc2019_Regev_samp5", "Bmsc2019_Regev_samp6", 
            "Bmsc2019_Regev_b1", "Bmsc2019_Regev_b2", "Bmsc2019_Regev_b3", 
            "Bmsc2019_Regev_b4", "Bmsc2019_Regev_bm1", "Bmsc2019_Regev_bm2", 
            "Bmsc2019_Regev_bm3", "Bmsc2019_Regev_bm4", "BmscChondro_Long", 
            "BmscEndosteal_Ono_Fgfr3CE", "BmscEndosteal_Ono_Fgfr3CEp53cKO", 
            "BmscEndosteal_Ono_Fgfr3CEp53cHet", "BmscEndosteal_Ono_Gas1CE", 
            "BmscEndosteal_Ono_Prrx1cre18M", "BmscEndosteal_Ono_Prrx1creP21", 
            "BmscEndosteal_Ono_Prrx1creP21multiome", "BmscSpecification_Kishor_1", 
            "BmscSpecification_Kishor_2", "BmscSpecification_Kishor_3", "BmscSpecification_Kishor_4", 
            "BmscTime_Zhong_1M", "BmscTime_Zhong_1.5M", "BmscTime_Zhong_3M", 
            "BmscTime_Zhong_16M", "GrowthplateSox9_Abdul_WTP13", "GrowthplateSox9_Abdul_WTP19", 
            "PerichondrialP21_Matsushita_FR3CreCxcl12GfpP21", "PerichondrialP21_Matsushita_HesCreCxcl12GfpP21", 
            "Septoclasts_Kishor_Longbone1", "Septoclasts_Kishor_Longbone2", 
            "Septoclasts_Kishor_Pdgfra", "GrowthPlateP0_Mizuhashi", "Metaphysis_Yang_1", 
            "Metaphysis_Yang_2", "Metaphysis_Yang_3", "Metaphysis_Yang_4", 
            "Metaphysis_Yang_sorted")
merge1 <- mergeAll[,mergeAll$orig.ident%in%c(origin1)]
merge2 <- mergeAll[,!mergeAll$orig.ident%in%c(origin1)]
merge2@assays$originalexp <- NULL
merge1@assays$RNA <- merge1@assays$originalexp
merge1@assays$originalexp <- NULL
merge1@active.assay <- "RNA"
merge2@active.assay <- "RNA"

merge1@assays$RNA@key <- "RNA_"
merge2@assays$RNA@key

rm(mergeAll)
mergeAll <- merge(merge1,merge2)

for (i in unique(mergeAll$Sample)){
  if (length(mergeAll[,mergeAll$orig.ident==i]@assays$RNA@counts@i)==0){
    print(i)
  }
}
mergeAll <- mergeAll[,!mergeAll$Sample%in%c("SutureP0_Tower_WT","SutureP0_Tower_Mut")]
SaveH5Seurat(mergeAll,"../unimportant_processed_data/8.25_integrateAll/8.26_merge.h5Seurat",overwrite = T)
Convert("../unimportant_processed_data/8.25_integrateAll/8.26_merge.h5Seurat",dest="h5ad",,overwrite = T)
