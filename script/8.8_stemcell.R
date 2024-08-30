stemCells <- readRDS("3.9_wt_integrate//processed_data/4.22_highlightcell.Rds")
stemCells <- c(stemCells,list(`Sox9+ chondrocyte` = sc_sox9_cell, `Hoxa11 Mesenchyme` = sc_hoxa11_cell))
stemCellsSelect <- stemCells[c("Acta2+ BMSC",  "Axin2+ MSC", 
                         "Clec11a+ BMSC", "Hypertrophic chondrocyte", "Chondrocyte", 
                         "Early perichondrial cell", "Fox2a+ chondrocyte", "Gli1+ MSC", 
                         "Grem1+ BMSC", "Hhip+ MSC", "Cd168+ skeletal stem/progenitor cell", 
                         "Lepr+ BMSC", "Mcam+ BMSC", "Msx2+ MSC", "Pdgfra+ BMSC", "Prrx1+ MSC", 
                         "Pthlh+ chondrocyte", "pvSSC", "ocSSC", "Hes1+ Perichondrium","Sox9+ chondrocyte","Hoxa11 Mesenchyme"
)]

stemCellsLow <- stemCells["Cd168+ skeletal stem/progenitor cell"]
stemCellsMiddle <- stemCells[c("pvSSC", "ocSSC","Mcam+ BMSC","Pdgfra+ BMSC")]
stemCellsHigh <- setdiff(stemCellsSelect,stemCellsMiddle)
stemCellsHigh <-  setdiff(stemCellsHigh,stemCellsLow)

stemCellsMiddle <- stemCellsMiddle %>% unlist() %>% unique
stemCellsHigh <- stemCellsHigh %>% unlist() %>% unique
stemCellsLow <- stemCellsLow %>% unlist() %>% unique
myCol <- RColorBrewer::brewer.pal(3,"Set1")
p1=DimPlot(wt_integrate,cells.highlight = stemCellsHigh,cols.highlight = myCol[1])
p2=DimPlot(wt_integrate,cells.highlight = stemCellsMiddle,cols.highlight =  myCol[2])
p3=DimPlot(wt_integrate,cells.highlight = stemCellsLow,cols.highlight =  myCol[3])
p1|p2|p3
dir.create("revision/results/opc")
ggsave("revision/results/opc/confident_level.pdf",width = 12,height = 3)
saveRDS(stemCells,"3.9_wt_integrate//processed_data/8.8_highlightcell.Rds")
