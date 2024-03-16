#== benchmark plot----------------------------
#== plot gradient------------------
lamian_composition_p
traj_grad <- read.csv("processed_data/11.23_process_null/11.23_trajdiff_gradiant_metrics.csv",row.names = 1)
condiment_grad <- read.csv("processed_data/11.23_process_null/11.23_condiment_metrics.csv",row.names = 1)
p_val_list <- data.frame("lamian"=lamian_composition_p$p_val,"trajDiff"=traj_grad$p,"condiments"=condiment_grad$p_val)
p_val_list$fold <- seq(0,1,0.1)
long_p <- gather(p_val_list, key = "methods", value = "p_val",-fold)
long_p$fold=1-long_p$fold
long_p$log_p=-log(long_p$p_val+0.01)
ggplot(data = long_p,aes(x =fold , y=log_p,color=methods))+
  geom_point()+
  geom_line()+theme_bw()
dir.create("result/11.24_benchmark_composition")
ggsave("result//11.24_benchmark_composition/11.24_grad_p.pdf",width = 6,height = 6)
write.csv(long_p,"processed_data/11.23_process_null/11.23_long_p_grad.csv")
stat_list <- data.frame("lamian"=lamian_composition_p$z_val,"trajDiff"=traj_grad$ratio,"condiments"=condiment_grad$statics)
min_max <- function(x){
  (x-min(x))/(max(x)-min(x))
}
stat_list <- sapply(stat_list,min_max)%>%as.data.frame()
stat_list$fold <- seq(0,1,0.1)
long_stat <- gather(stat_list, key = "methods", value = "stat",-fold)
long_stat$fold=1-long_stat$fold
ggplot(data = long_stat,aes(x =fold, y=stat,color=methods))+
  geom_point()+
  geom_line()+theme_bw()
dir.create("result/11.24_benchmark_composition")
ggsave("result//11.24_benchmark_composition/11.24_grad_stat.pdf",width = 6,height = 6)


#== plot null-------------------------
traj_null <- read.csv("processed_data/11.23_process_null/11.24_traj_null.csv",row.names = 1)
lamian_composition_null
condiment_metrics_shuffle

p_val_list <- data.frame("lamian"=lamian_composition_null$p_val,"trajDiff"=traj_null$p,"condiments"=condiment_metrics_shuffle2$p_val)
p_val_list$trial <- seq(1,10,1)
long_p <- gather(p_val_list, key = "methods", value = "p_val",-trial)

long_p$log_p=-log(long_p$p_val+0.01)
ggplot(data = long_p,aes(x =trial , y=log_p,color=methods))+
  geom_point()+
  geom_line()+theme_bw()

ggsave("result//11.24_benchmark_composition/11.24_null_p.pdf",width = 6,height = 6)


score_1_1 <- 1-sum(traj_null$p<0.05)/10
score_1_2 <- 1-sum(lamian_composition_null$p<0.05)/10
score_1_3 <- 1-sum(condiment_metrics_shuffle2$p_val<0.05)/10
score_2_1 <- sum(traj_grad$p<0.05)/11
score_2_2 <- sum(lamian_composition_p$p<0.05)/11
score_2_3 <- sum(condiment_grad$p<0.05)/11
score_df <- data.frame("null" =c(score_1_1,score_1_2,score_1_3),
                       "acc" = c(score_2_1,score_2_2,score_2_3))
rownames(score_df) <- c("trajDiff","lamian","condiments")
score_df$methods <- rownames(score_df)
long_score <- gather(score_df, key = "type", value = "score",-methods)
  
ggplot(data = score_df,aes(x =null , y=acc,color=methods))+ylim(0,1)+xlim(0,1)+
  geom_point(size=5)+  geom_text(aes(label = methods), nudge_y = 0.05, check_overlap = TRUE) +
  theme_bw()
ggsave("result//11.24_benchmark_composition/11.24_metric_score.pdf",width = 8,height = 6)
write.csv(score_df,"processed_data/11.23_process_null/11.23_score_df.csv")
write.csv(p_val_list,"processed_data/11.23_process_null/11.23_p_val_null.csv")
write.csv(p_val_list,"processed_data/11.23_process_null/11.23_p_val_null2.csv")
