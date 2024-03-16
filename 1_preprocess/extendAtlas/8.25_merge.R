rm(list=ls())
library(Seurat)
library(Matrix)
setwd("../../../disk1/limb/8.25_full_integration/")
wtintegrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
wtintegrate <- DietSeurat(wtintegrate,data = FALSE,scale.data = FALSE)
empty_matrix <- sparseMatrix(dims = c(nrow(wtintegrate),ncol(wtintegrate)), i={}, j={})
empty_matrix <- as(empty_matrix, "dgCMatrix")
dimnames(empty_matrix) <- dimnames(wtintegrate)
wtintegrate <- SetAssayData(wtintegrate, slot = "data", new.data = empty_matrix)
wtintegrate@assays$originalexp@data <- NULL

unique(wtintegrate$Project)

metadata <- readxl::read_xlsx("data/8.18_sampleMeta.xlsx")
project <- unique(metadata$Project)
setdiff(project,unique(wtintegrate$Project))

#= LimbFgf23_Ayako-----------
LimbFgf23_Ayako <- readRDS("../data/cleandata/LimbFgf23_Ayako.Rds")
LimbFgf23_Ayako <- DietSeurat(LimbFgf23_Ayako)

#== Articular_Arce------------------
Articular_Arce <- readRDS("../data/cleandata/Articular_Arce.Rds")
Articular_Arce <- DietSeurat(Articular_Arce)
unique(Articular_Arce$orig.ident)

#== Articular_Sebastian--------------------
Articular_Sebastian <- readRDS("../data/cleandata/Articular_Sebastian.Rds")
Articular_Sebastian <- DietSeurat(Articular_Sebastian)
unique(Articular_Sebastian$orig.ident)

#== BMSC-Stress_Mo----------------------
BMSC_Stress_Mo <- readRDS("../data/cleandata/BMSC-Stress_Mo.Rds")
BMSC_Stress_Mo <- DietSeurat(BMSC_Stress_Mo)
unique(BMSC_Stress_Mo$orig.ident)
new.id <- c("BMSC-Stress_Mo_Adult1",
            "BMSC-Stress_Mo_Adult2",
            "BMSC-Stress_Mo_Aging",
            "BMSC-Stress_Mo_Fracture",
            "BMSC-Stress_Mo_Irradiation",
            "BMSC-Stress_Mo_PrxBM",
            "BMSC-Stress_Mo_Rosiglitazone"
)
names(new.id) <- unique(BMSC_Stress_Mo$orig.ident)
Idents(BMSC_Stress_Mo) <- BMSC_Stress_Mo$orig.ident
BMSC_Stress_Mo <- RenameIdents(BMSC_Stress_Mo,new.id)
BMSC_Stress_Mo$orig.ident <- Idents(BMSC_Stress_Mo)

#== BmscAging_Young------------------
BmscAging_Young <- readRDS("../data/cleandata/BmscAging_Young.Rds")
BmscAging_Young <- DietSeurat(BmscAging_Young)
unique(BmscAging_Young$orig.ident)

#==BmscDiffereation_Wolock----------------------
BmscDiffereation_Wolock <- readRDS("../data/cleandata/BmscDiffereation_Wolock")
BmscDiffereation_Wolock <- DietSeurat(BmscDiffereation_Wolock)
unique(BmscDiffereation_Wolock$orig.ident)
BmscDiffereation_Wolock$orig.ident <- "BmscDiffereation_Wolock"

#==BmscDll4_Xu-------------------------
BmscDll4_Xu <- readRDS("../data/cleandata/BmscDll4_Xu.Rds")
BmscDll4_Xu <- DietSeurat(BmscDll4_Xu)
unique(BmscDll4_Xu$orig.ident)

#==BmscEpi_Kanazawa------------------------
BmscEpi_Kanazawa <- readRDS("../data/cleandata/BmscEpi_Kanazawa.Rds")
BmscEpi_Kanazawa <- DietSeurat(BmscEpi_Kanazawa)
unique(BmscEpi_Kanazawa$orig.ident)
BmscEpi_Kanazawa$orig.ident <- "BmscEpi_Kanazawa"

#==BmscFac2019_Anastasia--------------------
BmscFac2019_Anastasia <- readRDS("../data/finalData//BmscFac2019_Anastasia.Rds")
BmscFac2019_Anastasia <- DietSeurat(BmscFac2019_Anastasia)
unique(BmscFac2019_Anastasia$orig.ident)
new.id <- c("BmscFac2019_Anastasia_5FUcol23",
            "BmscFac2019_Anastasia_5FUlepr",
            "BmscFac2019_Anastasia_5FUvecad",
            "BmscFac2019_Anastasia_CNTRL",
            "BmscFac2019_Anastasia_CNTRLcol23",
            "BmscFac2019_Anastasia_NICHEcol23",
            "BmscFac2019_Anastasia_NICHElepr",
            "BmscFac2019_Anastasia_NICHEvecad")
names(new.id) <-unique(BmscFac2019_Anastasia$orig.ident)
Idents(BmscFac2019_Anastasia) <- BmscFac2019_Anastasia$orig.ident
BmscFac2019_Anastasia <- RenameIdents(BmscFac2019_Anastasia,new.id)
BmscFac2019_Anastasia$orig.ident <- Idents(BmscFac2019_Anastasia)

#==BmscGli_Shi---------------------------
BmscGli_Shi <- readRDS("../data/cleandata/BmscGli_Shi.Rds")
BmscGli_Shi <- DietSeurat(BmscGli_Shi)
BmscGli_Shi$orig.ident <- "BmscGli_Shi"

#== BmscGVHD_Gu--------------
BmscGVHD_Gu <- readRDS("../data/cleandata/BmscGVHD_Gu.Rds")
BmscGVHD_Gu <- DietSeurat(BmscGVHD_Gu)
unique(BmscGVHD_Gu$orig.ident)

#== BmscHuman_Jardine--------------
BmscHuman_Jardine <- readRDS("../data/cleandata/BmscHuman_Jardine.Rds")
BmscHuman_Jardine <- DietSeurat(BmscHuman_Jardine)
unique(BmscHuman_Jardine$orig.ident)

#== BmscInjury_Zhong-----------------
BmscInjury_Zhong <- readRDS("../data/finalData/BmscInjury_Zhong.Rds")
BmscInjury_Zhong <- DietSeurat(BmscInjury_Zhong)
unique(BmscInjury_Zhong$orig.ident)
BmscInjury_Zhong$orig.ident <- "BmscInjury_Zhong"

#== BmscMds_Hayashi-------------------
BmscMds_Hayashi <- readRDS("../data/cleandata/BmscMds_Hayashi.Rds")
BmscMds_Hayashi <- DietSeurat(BmscMds_Hayashi)
unique(BmscMds_Hayashi$orig.ident)


#== BmscMyelofibrosis_Leimkühler-------------------
BmscMyelofibrosis_Leimkühler <- readRDS("../data/cleandata/BmscMyelofibrosis_Leimkühler.Rds")
BmscMyelofibrosis_Leimkühler <- DietSeurat(BmscMyelofibrosis_Leimkühler)
unique(BmscMyelofibrosis_Leimkühler$orig.ident)

#== BMSCNiche2019_Baccin-------------------
BMSCNiche2019_Baccin <- readRDS("../data/finalData/BMSCNiche2019_Baccin.Rds")
BMSCNiche2019_Baccin <- DietSeurat(BMSCNiche2019_Baccin)
BMSCNiche2019_Baccin$orig.ident <- "BMSCNiche2019_Baccin"
unique(BMSCNiche2019_Baccin$orig.ident)

#== BmscPth1r_Fu--------------------------
BmscPth1r_Fu <- readRDS("../data/cleandata/BmscPth1r_Fu.Rds")
BmscPth1r_Fu <- DietSeurat(BmscPth1r_Fu)
unique(BmscPth1r_Fu$orig.ident)

#== BmscSp7Cre_Youlten------------------------
BmscSp7Cre_Youlten <- readRDS("../data/cleandata/BmscSp7KO_Youlten.Rds")
BmscSp7Cre_Youlten <- DietSeurat(BmscSp7Cre_Youlten)
unique(BmscSp7Cre_Youlten$orig.ident)

#===CalvariaHuman_He-----------------
CalvariaHuman_He_1 <- readRDS("../data/finalData//CalvariaHuman_He_1.Rds")
CalvariaHuman_He_2 <- readRDS("../data/finalData//CalvariaHuman_He_2.Rds")
CalvariaHuman_He_1 <- DietSeurat(CalvariaHuman_He_1)
CalvariaHuman_He_1$orig.ident <- "CalvariaHuman_He_1"
CalvariaHuman_He_2 <- DietSeurat(CalvariaHuman_He_2)
CalvariaHuman_He_2$orig.ident <- "CalvariaHuman_He_2"
unique(BmscSp7Cre_Youlten$orig.ident)

#== CranioRepairP75_Xu-------------
CranioRepairP75_Xu <- readRDS("../data/cleandata/CranioRepairP75_Xu.Rds")
CranioRepairP75_Xu <- DietSeurat(CranioRepairP75_Xu)

#== Digit2019_Storer--------------------
Digit2019_Storer_E11 <- readRDS("../data/finalData//Digit2019_Storer_E11.Rds")
Digit2019_Storer_E14 <- readRDS("../data/finalData//Digit2019_Storer_E14.Rds")
Digit2019_Storer_Gen <- readRDS("../data/finalData//Digit2019_Storer_Gen.Rds")
Digit2019_Storer_NonGen <- readRDS("../data/finalData//Digit2019_Storer_NonGen.Rds")
Digit2019_Storer_P3 <- readRDS("../data/finalData//Digit2019_Storer_P3.Rds")
Digit2019_Storer_E11 <- DietSeurat(Digit2019_Storer_E11)
Digit2019_Storer_E14 <- DietSeurat(Digit2019_Storer_E14)
Digit2019_Storer_NonGen <- DietSeurat(Digit2019_Storer_NonGen)
Digit2019_Storer_Gen <- DietSeurat(Digit2019_Storer_Gen)
Digit2019_Storer_P3 <- DietSeurat(Digit2019_Storer_P3)
Digit2019_Storer_E11$orig.ident <- "Digit2019_Storer_E11"
Digit2019_Storer_E14$orig.ident <- "Digit2019_Storer_E14"
Digit2019_Storer_Gen$orig.ident <- "Digit2019_Storer_Gen"
Digit2019_Storer_NonGen$orig.ident <- "Digit2019_Storer_NonGen"
Digit2019_Storer_P3$orig.ident <- "Digit2019_Storer_P3"

new.id <- c(
            "Digit2019_Storer_Gen7DPA",
            "Digit2019_Storer_Gen14DPAcreDmp1",
            "Digit2019_Storer_UninjuredcreDmp1",
            "Digit2019_Storer_Gen10DPA",
            "Digit2019_Storer_Gen14DPA1",
            "Digit2019_Storer_Gen14DPA2",
            "Digit2019_Storer_Gen28DPA",
            "Digit2019_Storer_Gen56DPA",
            "Digit2019_Storer_Uninjured1",
            "Digit2019_Storer_Uninjured2"
)
names(new.id) <- unique(Digit2019_Storer_Gen$orig.ident)
Idents(Digit2019_Storer_Gen) <- Digit2019_Storer_Gen$orig.ident
Digit2019_Storer_Gen <- RenameIdents(Digit2019_Storer_Gen,new.id)
Digit2019_Storer_Gen$orig.ident <- Idents(Digit2019_Storer_Gen)

new.id <- c("Digit2019_Storer_NonGen10DPA2",
            "Digit2019_Storer_NonGen14DPA",
            "Digit2019_Storer_NonGen10DPA1"
)

names(new.id) <- unique(Digit2019_Storer_NonGen$orig.ident)
Idents(Digit2019_Storer_NonGen) <- Digit2019_Storer_NonGen$orig.ident
Digit2019_Storer_NonGen <- RenameIdents(Digit2019_Storer_NonGen,new.id)
Digit2019_Storer_NonGen$orig.ident <- Idents(Digit2019_Storer_NonGen)


#== DigitLmx1b_Miller-----------------------
DigitLmx1b_Miller <- readRDS("../data/cleandata/DigitLmx1b_Miller.Rds")
DigitLmx1b_Miller <- DietSeurat(DigitLmx1b_Miller)

mergeSeurat <- merge(x = wtintegrate, y = list(Articular_Arce, Articular_Sebastian, BMSC_Stress_Mo,
                                               BmscAging_Young, BmscDiffereation_Wolock, BmscDll4_Xu,
                                               BmscEpi_Kanazawa, BmscFac2019_Anastasia, BmscGli_Shi, BmscGVHD_Gu,
                                               BmscHuman_Jardine, BmscInjury_Zhong, BmscMds_Hayashi, BmscMyelofibrosis_Leimkühler,
                                               BMSCNiche2019_Baccin, BmscPth1r_Fu, BmscSp7Cre_Youlten,
                                               CalvariaHuman_He_1, CalvariaHuman_He_2, CranioRepairP75_Xu,
                                               Digit2019_Storer_E11, Digit2019_Storer_E14, Digit2019_Storer_Gen,
                                               Digit2019_Storer_NonGen, Digit2019_Storer_P3, DigitLmx1b_Miller,LimbFgf23_Ayako))
dput(ls())


saveRDS(mergeSeurat,"../unimportant_processed_data/8.25_integrate1.Rds")
setdiff(project,unique(mergeSeurat$Project))

#== DigitRegeneration_Johnson--------------
DigitRegeneration_Johnson <- readRDS("../data/cleandata/DigitRegeneration_Johnson.Rds")
DigitRegeneration_Johnson <- DietSeurat(DigitRegeneration_Johnson)
unique(DigitRegeneration_Johnson$orig.ident)
new.id <- c("DigitRegeneration_Johnson_11dpa",
            "DigitRegeneration_Johnson_12dpa",
            "DigitRegeneration_Johnson_14dpa",
            "DigitRegeneration_Johnson_17dpa",
            "DigitRegeneration_Johnson_UA"
)
names(new.id) <-  unique(DigitRegeneration_Johnson$orig.ident)
Idents(DigitRegeneration_Johnson) <- DigitRegeneration_Johnson$orig.ident
DigitRegeneration_Johnson <- RenameIdents(DigitRegeneration_Johnson,new.id)
DigitRegeneration_Johnson$orig.ident <- Idents(DigitRegeneration_Johnson)

#== EnthesisGli_Fang------------------------
EnthesisGli_Fang <- readRDS("../data/cleandata/EnthesisGli_Fang.Rds")
EnthesisGli_Fang <- DietSeurat(EnthesisGli_Fang)
unique(EnthesisGli_Fang$orig.ident)
newid <- c("EnthesisGli_Fang_P11",
           "EnthesisGli_Fang_P18",
           "EnthesisGli_Fang_P56",
           "EnthesisGli_Fang_GliP7",
           "EnthesisGli_Fang_GliP11",
           "EnthesisGli_Fang_GliP14",
           "EnthesisGli_Fang_GliP18",
           "EnthesisGli_Fang_GliP21",
           "EnthesisGli_Fang_GliP28",
           "EnthesisGli_Fang_GliP42"
)
names(newid) <- unique(EnthesisGli_Fang$orig.ident)
Idents(EnthesisGli_Fang) <- EnthesisGli_Fang$orig.ident
EnthesisGli_Fang <- RenameIdents(EnthesisGli_Fang,newid)
EnthesisGli_Fang$orig.ident <- Idents(EnthesisGli_Fang)
#==Growthplate_Li-------------
Growthplate_Li <- readRDS("../data/finalData//Growthplate_Li.Rds")
Growthplate_Li <- DietSeurat(Growthplate_Li)
Growthplate_Li$orig.ident <- "Growthplate_Li"

#== HeterotopicActin_Mundy--------------------
HeterotopicActin_Mundy <- readRDS("../data/cleandata/HeterotopicActin_Mundy.Rds")
HeterotopicActin_Mundy <- DietSeurat(HeterotopicActin_Mundy)

#== LimbHumanDevelop-------------------------
LimbHumanDevelop_5w <- readRDS("../data/finalData/LimbHumanDevelop_He_5w.Rds")
LimbHumanDevelop_5w <- DietSeurat(LimbHumanDevelop_5w)
new.id <- c("LimbHumanDevelop_He_5w1",
            "LimbHumanDevelop_He_5w2",
            "LimbHumanDevelop_He_5w3"
)
names(new.id) <- unique(LimbHumanDevelop_5w$orig.ident)
Idents(LimbHumanDevelop_5w) <- LimbHumanDevelop_5w$orig.ident
LimbHumanDevelop_5w <- RenameIdents(LimbHumanDevelop_5w,new.id)
LimbHumanDevelop_5w$orig.ident <- Idents(LimbHumanDevelop_5w)

LimbHumanDevelop_8w <- readRDS("../data/finalData/LimbHumanDevelop_He_8w.Rds")
LimbHumanDevelop_8w <- DietSeurat(LimbHumanDevelop_8w)
new.id <- c("LimbHumanDevelop_He_8w1",
            "LimbHumanDevelop_He_8w2",
            "LimbHumanDevelop_He_8w3"
)
names(new.id) <- unique(LimbHumanDevelop_8w$orig.ident)
Idents(LimbHumanDevelop_8w) <- LimbHumanDevelop_8w$orig.ident
LimbHumanDevelop_8w <- RenameIdents(LimbHumanDevelop_8w,new.id)
LimbHumanDevelop_8w$orig.ident <- Idents(LimbHumanDevelop_8w)

#== LimbImplant_VesPrey----------------
LimbImplant_VesPrey <- readRDS("../data/cleandata/_List.Rds")
LimbImplant_VesPrey_Acta2 <- LimbImplant_VesPrey$compliant
LimbImplant_VesPrey_Tmem100 <- LimbImplant_VesPrey$tmem
LimbImplant_VesPrey_Acta2$orig.ident <- "LimbImplant_VesPrey_Acta2"
LimbImplant_VesPrey_Tmem100$orig.ident <- "LimbImplant_VesPrey_Tmem100"
rm(LimbImplant_VesPrey)
LimbImplant_VesPrey_Acta2 <- DietSeurat(LimbImplant_VesPrey_Acta2)
LimbImplant_VesPrey_Tmem100 <- DietSeurat(LimbImplant_VesPrey_Tmem100)

#== LimbSSC_Chan----------------------
LimbSSC_Chan <- readRDS("../data/cleandata//LimbSSC_Chan.Rds")
LimbSSC_Chan <- DietSeurat(LimbSSC_Chan)
LimbSSC_Chan$orig.ident <- "LimbSSC_Chan"

#== MscMineCulture_Basel---------------------
MscMineCulture_Basel <- readRDS("../data/cleandata/MscMineCulture_Basel.Rds")
MscMineCulture_Basel <- DietSeurat(MscMineCulture_Basel)
unique(MscMineCulture_Basel$orig.ident)
MscMineCulture_Basel$orig.ident <- "MscMineCulture_Basel"

#== Periodontium_Nagata----------------------
Periodontium_Nagata <- readRDS("../data/annodata///Periodontium_Nagata.Rds")
Periodontium_Nagata <- DietSeurat(Periodontium_Nagata)

#==Periosteal2018_Shawon-----------------------
Periosteal2018_Shawon <- readRDS("../data/finalData/Periosteal2018_Shawon.Rds")
Periosteal2018_Shawon <- DietSeurat(Periosteal2018_Shawon)
Periosteal2018_Shawon$orig.ident <- "Periosteal2018_Shawon"

#== RibRegeneraton_Serowoky-----------------------
RibRegeneraton_Serowoky <- readRDS("../data/cleandata/RibRegeneraton_Serowoky.Rds")
RibRegeneraton_Serowoky <- DietSeurat(RibRegeneraton_Serowoky)

#== SkeletalMuscle_Julien---------------
SkeletalMuscle_Julien <- readRDS("../data/cleandata/SkeletalMuscle_Julien.Rds")
SkeletalMuscle_Julien <- DietSeurat(SkeletalMuscle_Julien)


#==Sp7Cre2019_Bohm--------------------
Sp7Cre2019_Bohm <- readRDS("../data/finalData//Sp7Cre2019_Bohm.Rds")
Sp7Cre2019_Bohm <- DietSeurat(Sp7Cre2019_Bohm)
Sp7Cre2019_Bohm$orig.ident <- "Sp7Cre2019_Bohm"

#== SSC2021_Ambrosi-------------------------------
SSC2021_Ambrosi <- readRDS("../data/finalData/SSC2021_Ambrosi.Rds")
SSC2021_Ambrosi <- DietSeurat(SSC2021_Ambrosi)
SSC2021_Ambrosi$orig.ident <- "SSC2021_Ambrosi"

#== SutureP0_Tower--------------------
SutureP0_Tower_WT <- readRDS("../data/cleandata/SutureP0_Tower_WT.Rds")
SutureP0_Tower_Mut <- readRDS("../data/cleandata/SutureP0_Tower_Mut.Rds")
SutureP0_Tower_Mut <- DietSeurat(SutureP0_Tower_Mut,assays = )
SutureP0_Tower_Mut@images <- list()
SutureP0_Tower_Mut$orig.ident <- "SutureP0_Tower_Mut"
SutureP0_Tower_WT$orig.ident <- "SutureP0_Tower_WT"
SutureP0_Tower_WT <- DietSeurat(SutureP0_Tower_WT)
SutureP0_Tower_WT@images <- list()

#== Mmp14_Chu-------------------------
Mmp14_Chu <- readRDS("../data/annodata//Mmp14_Chu.Rds")
Mmp14_Chu <- DietSeurat(Mmp14_Chu)

#== LimbBgn_Shainer-------------------
LimbBgn_Shainer <- readRDS("../data/annodata/LimbBgn_Shainer.Rds")
LimbBgn_Shainer <- DietSeurat(LimbBgn_Shainer)
setdiff(project,unique(mergeSeurat$orig.ident))
setdiff(metadata$Sample,unique(mergeSeurat$orig.ident))

#== Frontal_Holmes--------------
Frontal_Holmes <- readRDS("../data/annodata/Frontal_Holmes_KO.Rds")
Frontal_Holmes <- DietSeurat(Frontal_Holmes)

#==Bmsc2019_Regev_Leu------------------
Bmsc2019_Regev_Leu <- readRDS("../data/cleandata/Bmsc2019_Regev_Leu.Rds")
new.id <- c(
            "Bmsc2019_Regev_samp1",
            "Bmsc2019_Regev_samp2",
            "Bmsc2019_Regev_samp3",
            "Bmsc2019_Regev_samp4",
            "Bmsc2019_Regev_samp5",
            "Bmsc2019_Regev_Leu1",
            "Bmsc2019_Regev_Leu2",
            "Bmsc2019_Regev_Leu3",
            "Bmsc2019_Regev_Leu4"
)
names(new.id) <- unique(Bmsc2019_Regev_Leu$orig.ident)
Idents(Bmsc2019_Regev_Leu) <- Bmsc2019_Regev_Leu$orig.ident
Bmsc2019_Regev_Leu <- RenameIdents(Bmsc2019_Regev_Leu,new.id)
Bmsc2019_Regev_Leu$orig.ident <- Idents(Bmsc2019_Regev_Leu)

Bmsc2019_Regev_Leu <- DietSeurat(Bmsc2019_Regev_Leu)

#==BmscSpecification_Kishor------------------
BmscSpecification_Kishor_Pdgfa <- readRDS("../data/annodata/BMSCSpecification_Kishor_Pdgfa.Rds")
BmscSpecification_Kishor_Pdgfa <- DietSeurat(BmscSpecification_Kishor_Pdgfa)


#==CranioSoxc_Angelozzi------------------
CranioSoxc_Angelozzi <- readRDS("../data/annodata//CranioSoxc_Angelozzi_KO.Rds")
CranioSoxc_Angelozzi <- DietSeurat(CranioSoxc_Angelozzi)
new.id <- c("CranioSoxc_Angelozzi_Prx1CreE11.5",
            "CranioSoxc_Angelozzi_Prx1CreE12.5",
            "CranioSoxc_Angelozzi_Prx1CreE13.5",
            "CranioSoxc_Angelozzi_Prx1CreE15.5",
            "CranioSoxc_Angelozzi_Prx1CreE17.5"
)

names(new.id) <- unique(CranioSoxc_Angelozzi$orig.ident)
Idents(CranioSoxc_Angelozzi) <- CranioSoxc_Angelozzi$orig.ident
CranioSoxc_Angelozzi <- RenameIdents(CranioSoxc_Angelozzi,new.id)
CranioSoxc_Angelozzi$orig.ident <- Idents(CranioSoxc_Angelozzi)
unique(mergeSeurat$orig.ident)


#== LimbG610C_Gorrell-------------
LimbG610C_Gorrell <- readRDS("../data/annodata/LimbG610C_Gorrell_KO.Rds")
LimbG610C_Gorrell <- DietSeurat(LimbG610C_Gorrell)
dim(LimbG610C_Gorrell)

#== GrowthplateSox9_Abdul_KO-------------------------
GrowthplateSox9_Abdul_KO <- readRDS("../data/annodata/GrowthplateSox9_Abdul_KO.Rds")
GrowthplateSox9_Abdul_KO <- DietSeurat(GrowthplateSox9_Abdul_KO)
unique(GrowthplateSox9_Abdul_KO$orig.ident)


mergeSeurat$orig.ident[mergeSeurat$orig.ident=="CB22"]="CalvariaHuman_He_1"

#== Forelimb_He_Smart--------------
Forelimb_He_Smart <- readRDS("../data/cleandata/Forelimb_He_smart.Rds")
Forelimb_He_Smart <- DietSeurat(Forelimb_He_Smart)
unique(Forelimb_He_Smart$orig.ident)

mergeSeurat$orig.ident[mergeSeurat$orig.ident=="Digit2019_Storer_Gen"] <- Digit2019_Storer_Gen$orig.ident
mergeSeurat$orig.ident[mergeSeurat$orig.ident=="Digit2019_Storer_NonGen"] <- Digit2019_Storer_NonGen$orig.ident


#== Msx------------------------
MsxCranio_Zhang <- readRDS("../preprocess/11.28_Msx/result/12.3_MsxSub.Rds")
DimPlot(MsxCranio_Zhang)
new.id <- c("MsxCranio_Zhang_W1Con1",
            "MsxCranio_Zhang_W1N1",
            "MsxCranio_Zhang_W2Con1",
            "MsxCranio_Zhang_W2Con2",
            "MsxCranio_Zhang_W2N1",
            "MsxCranio_Zhang_W2N2"
            )
names(new.id) <- unique(MsxCranio_Zhang$orig.ident)
Idents(MsxCranio_Zhang) <- MsxCranio_Zhang$orig.ident
MsxCranio_Zhang <- RenameIdents(MsxCranio_Zhang,new.id)
MsxCranio_Zhang$orig.ident <- Idents(MsxCranio_Zhang)
MsxCranio_Zhang <- DietSeurat(MsxCranio_Zhang)
merge2 <- merge(Bmsc2019_Regev_Leu,list(BmscSpecification_Kishor_Pdgfa, CranioSoxc_Angelozzi,
                  DigitRegeneration_Johnson, EnthesisGli_Fang, Forelimb_He_Smart,
                  Frontal_Holmes, Growthplate_Li, GrowthplateSox9_Abdul_KO,
                  HeterotopicActin_Mundy, LimbBgn_Shainer, LimbG610C_Gorrell,
                  LimbHumanDevelop_5w, LimbHumanDevelop_8w, LimbImplant_VesPrey_Acta2,
                  LimbImplant_VesPrey_Tmem100, LimbSSC_Chan, 
                  Mmp14_Chu, MscMineCulture_Basel, MsxCranio_Zhang,
                  Periodontium_Nagata, Periosteal2018_Shawon,
                  RibRegeneraton_Serowoky, SkeletalMuscle_Julien,
                  Sp7Cre2019_Bohm, SSC2021_Ambrosi, SutureP0_Tower_Mut, SutureP0_Tower_WT
))

saveRDS(merge2,"../unimportant_processed_data/8.25_integrate2.Rds")
rm(list=c("Bmsc2019_Regev_Leu","BmscSpecification_Kishor_Pdgfa", "CranioSoxc_Angelozzi",
           "DigitRegeneration_Johnson", "EnthesisGli_Fang", "Forelimb_He_Smart",
           "Frontal_Holmes", "Growthplate_Li", "GrowthplateSox9_Abdul_KO",
           "HeterotopicActin_Mundy", "LimbBgn_Shainer", "LimbG610C_Gorrell",
           "LimbHumanDevelop_5w", "LimbHumanDevelop_8w", "LimbImplant_VesPrey_Acta2",
           "LimbImplant_VesPrey_Tmem100", "LimbSSC_Chan", 
           "Mmp14_Chu", "MscMineCulture_Basel", "MsxCranio_Zhang",
           "Periodontium_Nagata", "Periosteal2018_Shawon",
           "RibRegeneraton_Serowoky", "SkeletalMuscle_Julien",
           "Sp7Cre2019_Bohm", "SSC2021_Ambrosi", "SutureP0_Tower_Mut", "SutureP0_Tower_WT"
))
mergeAll <- merge(mergeSeurat,merge2)
saveRDS(mergeAll,"../important_processed_data/8.26_merge_71w.Rds")
