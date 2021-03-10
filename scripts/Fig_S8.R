library(pheatmap)
library(ComplexHeatmap)
library("gplots")

hk_arg_id1_ACT_sample_new<-read.table('./HK_args_class_freq.txt',header = T)
ggplot(data = hk_arg_id1_ACT_sample_new)+
  geom_bar(mapping = aes(x=sample,fill=ACT),position = "fill")+
  theme(axis.text.x=element_text(angle = 90, size=12,hjust = 1),
        axis.text.y=element_text(size=12))
cols_column<-read.table('./HK_74_ACT_ARG_new_col.txt')
cols_column
cols_column<-cols_column[,1]
cols_column
level_order <- c(paste0('HK',1:74))
level_order
p1<-ggplot(data = hk_arg_id1_ACT_sample_new)+
  geom_bar(mapping = aes(x=factor(sample, level = level_order),fill=ACT),position = "fill")+
  theme(axis.text.x=element_text(angle = 90, size=12,hjust = 1,colour = cols_column),
        axis.text.y=element_text(size=12),axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),legend.text=element_text(size=10))+
  labs(x="Subjects from Hong Kong",y="Frequency of each class of ARGs")+
  scale_fill_manual(values = pal_simpsons('springfield')(16)[3:15])
p1

#### plot p2 tax abundance in HK
HK_merged_abundance_table_species<-read.table('./merged_abundance_table_species_table_new.txt',header = T, row.names = 1)
HK_merged_abundance_table_species_log10<-log10(HK_merged_abundance_table_species)
HK_samples<-paste0('HK',1:21)
HK_samples
cols_column <- rep('black', ncol(HK_merged_abundance_table_species_log10))
cols_column[colnames(HK_merged_abundance_table_species_log10) %in% HK_samples] <- 'red'
names(cols_column) <- colnames(HK_merged_abundance_table_species_log10)

cols_column <- rep('HK_Yu', ncol(HK_merged_abundance_table_species_log10))
cols_column[colnames(HK_merged_abundance_table_species_log10) %in% HK_samples] <- 'HK_HKU'
annotation_col = data.frame(Group=cols_column,row.names = colnames(HK_merged_abundance_table_species_log10))

mycolors <- pal_npg(c("nrc"))(10)[c(1,4)]
names(mycolors) <- unique(annotation_col$Group)
mycolors <- list(Group = mycolors)

HK_merged_abundance_table_species_log10[is.infinite(as.matrix(HK_merged_abundance_table_species_log10))]<-0

heatmap.2(as.matrix(HK_merged_abundance_table_species_log10), scale = "none",
          col =  col,trace = "none", density.info = "none",margins = c(5, 14),keysize = 1.1,colCol = cols_column,cexCol = 0.9)

row.names(HK_merged_abundance_table_species_log10) <- gsub("_"," ",row.names(HK_merged_abundance_table_species_log10))
p2<-pheatmap(HK_merged_abundance_table_species_log10,fontsize = 10,
             color=(colorRampPalette(brewer.pal(9, "RdYlBu"))(1000)),
             annotation_col = annotation_col,
             annotation_colors = mycolors,
             silent = TRUE)
p2

mycolors <- colorRampPalette(brewer.pal(9,'Blues'))(length(unique(annotation_col$Region)))
names(mycolors) <- unique(annotation_col$Region)
mycolors <- list(Region = mycolors)

## plot HK_IS_ARGs_ratio
hk_arg_with_IS_ratio<-read.table('./HK_IS_ARGs_ratio_new.txt',header = T)
hk_arg_with_IS_ratio$Age <- as.factor(hk_arg_with_IS_ratio$Age)
hk_arg_with_IS_ratio$sample = factor(hk_arg_with_IS_ratio$sample, 
                                     levels=hk_arg_with_IS_ratio[order(hk_arg_with_IS_ratio$ratio), "sample"])
col <- ifelse(hk_arg_with_IS_ratio$Category == 0, "red", "black")
p3<-ggplot(hk_arg_with_IS_ratio,aes(y=sample,x=ratio,col=Age))+
  geom_point(size=8,alpha=0.8)+
  theme_bw()+
  theme(axis.text.y=element_text(size=12,colour=col),
                   axis.text.x = element_text(size = 12),
                   axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))+
  labs(y="Human subjects from Hong Kong",x="Frequency of ARGs with IS")+
  scale_color_manual(values = pal_npg(c("nrc"))(10)[c(4,3,9,1)])
p3<-p3+scale_fill_discrete(name="Age",labels = c("<18","18 > and <40", ">40 and <50","50 > and <60"))
p3


#### ratio of coverage
hk_all_info_cov<-read.table('./hk_all_info_cov2.txt',header =T,sep='\t',stringsAsFactors = F)
p4<-ggplot(hk_all_info_cov, aes(Gene, y = value, color = variable))+
  geom_point(aes(y = RPKM, col = "RPKM",size=2))+
  geom_point(aes(y = Coverage, col = "Coverage",size=2))+theme_bw()+
  theme(axis.text.x=element_text(),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text=element_text(size=12),legend.title=element_blank(),
        legend.key.size = unit(0.5, "cm"))+
  labs(x="ARGs with IS",y="Coverage of contigs/RPKM of ARGs")+
  scale_color_manual(values = pal_npg(c("nrc"))(10)[c(1,4)])
p4<-p4+guides(colour = guide_legend(show = FALSE))
p4


### hive for HK_ARGs

library(ape)
library(ggtree)
library(HiveR)
library(forcats)
#Merge data

HK_ARG_Tax_hive <- read.table("./HK_ARG_Tax_hive.txt",sep="\t",header = T)
HK_ARG_Tax_hive_node <- read.table("./HK_ARG_Tax_hive_node.txt",sep="\t",header = T)
hive_hk <- edge2HPD(HK_ARG_Tax_hive)
hive_hk$nodes$axis<-HK_ARG_Tax_hive_node$class
hive_hk$nodes[hive_hk$nodes$axis==1, ] <- arrange(hive_hk$nodes[hive_hk$nodes$axis==1, ], lab)
#hive_hk$edges$weight <- hive_hk$edges$weight*3-1
hive_hk$edges$weight<-log10(hive_hk$edges$weight)
hive_hk$nodes$size=2
hive_hk$edges$color <- ifelse(hive_hk$edges$weight<2, '#66ccff55','#ff990055')
tmp <- data.frame(node.lab=hive_hk$nodes$lab, node.text=hive_hk$nodes$lab, angle=0, radius=0, offset=0, hjust=1, vjust=0.5)
tmp$lab<-HK_ARG_Tax_hive_node$tag
hive_hk$nodes$color[c(5,9,12,15,16,17,19,20,21)]<-colors[13]
hive_hk$nodes$color[c(22,23,24)]<-colors[15]
hive_hk$nodes[hive_hk$nodes$axis=='1',]$color<-sapply(c(pal_simpsons('springfield')(14)), adjustcolor, alpha.f=0.9)

mutate(tmp,offset=ifelse(lab=='Tax'| lab=='ARG', -0.05, -0.03)) %>% ###
  select(-lab) %>%
  write.table('tmp_hive.txt', sep=',', quote=T, row.names = F, col.names = T)

p5<-plotHive(hive_hk, ch=0.2, method='ranknorm', bkgnd='white', anNodes = 'tmp_hive.txt',
             anNode.gpar=gpar(cex=1))


# cowplot::plot_grid(p1,p2$gtable,p3,p4,p5,labels = c('a.', 'b.','c.','d.','e'),label_size = 20,nrow = 3)
# left_col<- plot_grid(p1,p2$gtable,p3, labels = c('a.', 'b.',"c."), label_size = 20,ncol = 1,rel_heights =c(8,12,8))
# left_col
# right_col<-plot_grid(p4,p5, labels = c('d.', 'e.'), label_size = 20,ncol = 1,nrow = 3)
# right_col
# plot_grid(left_col, right_col,
#           rel_widths = c(3,1.5))

left_col<-plot_grid(plot_grid(p1),
                    NULL,
                    plot_grid(p2$gtable),nrow = 3,rel_heights = c(4,0.3,4))
left_col
top_row <- plot_grid(plot_grid(left_col),
                     NULL,
                     plot_grid(p3),nrow = 1,rel_widths=c(15,0.5,5))
top_row
sec_row<- plot_grid(p4,p5,rel_widths=c(6,10))
plot_grid(plot_grid(top_row),
          NULL,
          plot_grid(sec_row), ncol = 1,rel_heights = c(12,0.6,8))

