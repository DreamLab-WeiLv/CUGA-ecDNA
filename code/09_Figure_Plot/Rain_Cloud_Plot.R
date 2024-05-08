p_s1 <- 
  plot_data %>%
  ggplot(aes(x = feature, y = Value, color = feature, fill = feature)) +  
  geom_boxplot(width = 0.1, fill  = "transparent",size  = 0.4,outlier.shape = NA) +  
  geom_half_violin(position = position_nudge(x=0.1,y=0),
                   alpha = 0.8, side = 'top')+  
  geom_half_point(range_scale = 0.5,shape=8,
                  side  = "l",alpha = 0.6, size  = 2) + 
  scale_fill_manual(values = c('Other fSCNA'='#619CFF','ecDNA'='#d33839')) +  
  scale_color_manual(values = c('Other fSCNA'='#619CFF','ecDNA'='#d33839'))+
  scale_y_sqrt() +  
  facet_grid(~Index,scales = 'free')+
  coord_flip(clip = F) +  
  labs(x=NULL,y= "Index") + 
  theme_classic() + 
  stat_compare_means(comparisons = list(c("ecDNA","Other fSCNA")),
                     label = "p.format",method = "wilcox.test")+
  theme( 
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.background= element_rect(fill = "white", color = "white"),  
    axis.ticks.y = element_blank(),  
    plot.margin= margin(t = 10, r = 10, b = 10, l = 10), 
    axis.text = element_text(size = 10,color="black"),  
    legend.position = "none")
p_s1
dev.off()
pdf('Figure.pdf',width = 12,height = 5)
p_s1
dev.off()
