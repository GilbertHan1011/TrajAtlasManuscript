TRAV1_index <- RAVindex[,143]
names(TRAV1_index) <- rownames(RAVindex)
TRAV1_index <- TRAV1_index[TRAV1_index>0]
TRAV1_index <- TRAV1_index[order(TRAV1_index,decreasing = T)[1:20]]
TRAV1_index <- as.data.frame(TRAV1_index)
pdf("result/24.2.22_fig5_supp/overview.pdf",width = 2,height = 4)
Heatmap(TRAV1_index,cluster_rows = F,col=colorRamp2(c(0, 0.002, 0.004), c("DeepSkyBlue3", "white", "red")),rect_gp = gpar(col = "black"),)
dev.off()

TRAV1_index <- TRAV1_index[c(1:20)]%>%as.data.frame()
TRAV1_index <- RAVindex[,143]
names(TRAV1_index) <- rownames(RAVindex)
TRAV1_index <- TRAV1_index[TRAV1_index>0]
TRAV1_index <- TRAV1_index[order(TRAV1_index,decreasing = T)]
pdf("result/24.2.22_fig5_supp/overview2.pdf",width = 2,height = 4)
Heatmap(TRAV1_index,cluster_rows = F,col=colorRamp2(c(0, 0.002, 0.004), c("DeepSkyBlue3", "white", "red")),show_row_names = F,)
dev.off()

