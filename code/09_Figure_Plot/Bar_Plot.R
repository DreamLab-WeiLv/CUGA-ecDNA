d <- read.table('test.txt',sep = '\t',header = T)
p1 <- ggplot(d, mapping = aes(x= ecDNA, fill = HighestLevel))+ 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = list('LEVEL_1' = '#559e52','LEVEL_2'='#98b478','LEVEL_3A'='#90539b','LEVEL_3B'='#b698c2','LEVEL_4'='#3b3e3a',
                                  'Oncogenic'='#d1d2d1','None'='#d1d2d1'))+
  theme_classic()+
  labs(y = 'Highest level of actionability',
       x = 'ecDNA Status')+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid = element_blank())+
  coord_flip()
pdf('Figure.pdf',width = 6,height = 4)
print(p1)
dev.off()
