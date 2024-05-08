d <- as.data.frame(read_xlsx('info.xlsx',sheet = "SMG in ecDNA+ freq"))
d <- d[,c(1,4,5,6)]
p <- ggplot(d, aes(x = Freq_in_noec, y = Freq_in_ec)) +
  geom_point(size=5) +
  geom_abline(linetype=2,slope = 1)+ 
  geom_text_repel(aes(label = Gene), hjust = 0, vjust = 0,fontface='italic',size=3) +  
  scale_y_continuous(limits = c(-9,55),breaks = c(0,10,20,30,40,50,60))+
  scale_x_continuous(limits = c(-9,55),breaks = c(0,10,20,30,40,50,60))+
  labs(
    title = "Scatter plot of Freq_in_noec vs Freq_in_ec",
    x = "Freq_in_noec",
    y = "Freq_in_ec") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pdf('Figure.pdf',width = 5,height = 5)
print(p)
dev.off()
