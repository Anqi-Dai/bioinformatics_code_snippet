# This script was used to make count data of two genes in barplot and show the three batches as grouped

library('tidyverse')

# the data is normalized 
dat <- read_rds('../../data/whole_eset_60samples.RDS')
gExpr <- dat[featureNames(dat) %in% c('ENSG00000149346','ENSG00000164362'),]
# full names of the cell lines
CL.fulln <- c(rep('HUO3N1',3),
              rep('hFOB1.19',3),
              rep('HOS',3),
              rep('MG63',3),
              rep('SJSA1',3),
              rep('CAL72',3),
              rep('CAL78',3),
              rep('G292',3),
              rep('NOS',3),
              rep('NY',3),
              rep('SAOS2',3),
              rep('U2OS',3),
              rep('HUO9',3))


pData(gExpr)$source <- as.character(pData(gExpr)$source)
pData(gExpr)$source[1:39] <- CL.fulln

# format the df with the gene counts and the phenotype data

source.level <- unique(pData(gExpr)$source)

exprDF <- data.frame(
  Sample = sampleNames(dat),
  SLX4IP = t(exprs(gExpr))[,1],
  TERT = t(exprs(gExpr))[,2],
  Status = pData(gExpr)$group, 
  Source = pData(gExpr)$source
) %>%
  mutate(Status = as.character(Status)) %>%
  mutate(Status = ifelse(Status == 'ALTM', 'ALT-', ifelse(Status == 'ALTP', 'ALT+', 'PDX'))) %>%
  mutate(Status = factor(Status), Sample = as.character(Sample)) %>%
  mutate(Source = factor(Source, levels = source.level))


# make the full cell line name as the sample name
Col <- c('#00468B', '#EC0000','#42B440')

exprDF <- exprDF %>%
  mutate(Batch = rep(seq(1,3), 20)) %>%
  mutate(SampleFull = paste(Source, Batch, sep = '-'))


# make the new x coord of the plot (with gaps!)
Coord <- sort(c(seq(1, 77,4), seq(2, 78,4), seq(3, 79,4)))

exprDF.new <- exprDF %>%
  mutate(x = Coord) %>%
  dplyr::select(x, SampleFull, Status, SLX4IP, TERT) %>%
  gather(key = 'Gene', value = 'Cts' ,SLX4IP:TERT ) %>%
  mutate(Gene = as.factor(Gene))


# Use the ggplot to plot

ggplot(exprDF.new, aes(x = x, y = Cts,  color = Status, fill = Status)) +
  geom_bar(stat="identity", width = 0.7) +  
  scale_x_continuous(breaks = exprDF.new$x,
                     labels = exprDF.new$SampleFull,
                     expand=c(0,0)) +
  labs(x = '',
       y = 'Relative expression') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=4),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        # strip is the facet grid title
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(color = "black", size = 0.5),
        legend.position="bottom") + 
  scale_fill_manual(values = Col) +
  scale_color_manual(values = Col)  +
  facet_grid(. ~ Gene) +
  ggsave('../figs/gene.exprs.slx4ip.tert.barplot.jpg', width = 10, height =5, dpi = 300)
