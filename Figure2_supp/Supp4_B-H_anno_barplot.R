setwd("Desktop/disk1/limb/3.9_wt_integrate/")
library(tidyverse)
library(ggstance)
wt_integrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
unique(wt_integrate$paper_label)

wt_integrate$short_id <- wt_integrate$Project%>%
  stringr::str_split("_")%>%
  map(.,function(x)
    paste(substr(x,1,2),collapse = ""))%>%
  unlist()

wt_integrate$paper_label_id <- paste0(wt_integrate$short_id,"_",wt_integrate$paper_label)

paperIdCount <- table(wt_integrate$paper_label_id)%>%as.data.frame()

naId <- grep("NA",paperIdCount$Var1)
paperIdCount <- paperIdCount[-naId,]

meta <- wt_integrate@meta.data
leprBMSCTable <- meta[meta$C19_named=="Lepr+ BMSC",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%top_n(5)%>%
  mutate(Var1 = reorder(Var1, percent))
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/LeprBMSC.pdf",width = 5,height = 3)


leprBMSCTable <- meta[meta$C19_named=="Pre-ob",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%
  mutate(Var1 = reorder(Var1, percent))
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/LeprBMSC.pdf",width = 5,height = 3)



leprBMSCTable <- meta[meta$C19_named=="Ly6a+ MSC",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%
  mutate(Var1 = reorder(Var1, percent))%>%top_n(5)
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/Ly6a+MSC.pdf",width = 5,height = 3)


leprBMSCTable <- meta[meta$C19_named=="Fibro",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%
  mutate(Var1 = reorder(Var1, percent))%>%top_n(5)
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/Fibroblast.pdf",width = 5,height = 3)

leprBMSCTable <- meta[meta$C7_named=="Pericyte",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%
  mutate(Var1 = reorder(Var1, percent))%>%top_n(5)
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/Pericyte.pdf",width = 5,height = 3)


leprBMSCTable <- meta[meta$C19_named=="Hmmr+ CPC",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%
  mutate(Var1 = reorder(Var1, percent))%>%top_n(5)
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/CPC.pdf",width = 5,height = 3)

leprBMSCTable <- meta[meta$C19_named=="Mature Chondro",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%
  mutate(Var1 = reorder(Var1, percent))%>%top_n(5)
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/MatureChondro.pdf",width = 5,height = 3)

leprBMSCTable <- meta[meta$C19_named=="Late.MSC",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%
  mutate(Var1 = reorder(Var1, percent))%>%top_n(5)
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/LateMSC.pdf",width = 5,height = 3)

leprBMSCTable <- meta[meta$C19_named=="CPC",]
leprBMSCTableCount <- table(leprBMSCTable$paper_label_id)%>%as.data.frame()%>%
  right_join(paperIdCount,by="Var1")%>%
  mutate(percent=Freq.x/Freq.y)%>%
  arrange(desc(percent))%>%
  mutate(Var1 = reorder(Var1, percent))%>%top_n(5)
ggplot(leprBMSCTableCount, aes(x = percent, y = Var1,fill=percent)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = Var1,y=Var1,x=rep(0,nrow(leprBMSCTableCount))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.2.11_Fig2_supp/CPC.pdf",width = 5,height = 3)
