library(reshape2)
library(tidyverse)
library(tidyr)
library(dplyr)

all_arg_res_rpkm <- read.csv("../../microbiome/data/resfinder_all/1174all_arg_res_rpkm.txt",sep="\t", 
                             header = T, row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
head(all_arg_res_rpkm)
dim(all_arg_res_rpkm)
all_arg_res_rpkm_t<-t(all_arg_res_rpkm)
head(all_arg_res_rpkm_t)
all_arg_res_rpkm_t<-data.frame(all_arg_res_rpkm_t,check.names = FALSE)
#disttransform(all_arg_res_rpkm, method="hellinger")
dim(all_arg_res_rpkm_t)
all_arg_res_rpkm_t<-rownames_to_column(all_arg_res_rpkm_t, var = "sample") #### tibble is needed

sample_region<-read.csv('1174_sample.txt',header = F,sep='\t',row.names = 1,check.names = FALSE)
sample_region<-rownames_to_column(sample_region, var = "sample")
colnames(sample_region)<-c("sample",'region','type',"nation")
head(sample_region)
dim(sample_region)
write.csv(sample_region,file = "sample_region.csv",row.names = TRUE, col.names = TRUE,quote=FALSE)

all_arg_res_rpkm_with_region<-merge(sample_region,all_arg_res_rpkm_t,'sample')
dim(all_arg_res_rpkm_with_region)
head(all_arg_res_rpkm_with_region)
all_arg_res_rpkm_with_region$sum<-rowSums(all_arg_res_rpkm_with_region[,5:659])
all_arg_res_rpkm_with_region$log2sum<-log2(all_arg_res_rpkm_with_region$sum)
levels(all_arg_res_rpkm_with_region$region)
dim(all_arg_res_rpkm_with_region)

all_arg_res_rpkm_human<-subset(all_arg_res_rpkm_with_region,Host == "Human")
all_arg_res_rpkm_human<-all_arg_res_rpkm_human[,-c(2,3,4)]
rownames(all_arg_res_rpkm_human)<-all_arg_res_rpkm_human[,1]
all_arg_res_rpkm_human<-all_arg_res_rpkm_human[,-1]
dim(all_arg_res_rpkm_human)
all_arg_res_rpkm_human_clean<-all_arg_res_rpkm_human[, colSums(all_arg_res_rpkm_human != 0) > 0]
dim(all_arg_res_rpkm_human_clean)
all_arg_human_df<-data.frame(colnames(all_arg_res_rpkm_human_clean))
colnames(all_arg_human_df)<-'arg'

arg_ACT<-read.csv('data/arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
library(tibble)
arg_ACT<-rownames_to_column(arg_ACT, var = "arg")
colnames(arg_ACT)<-c("arg",'ACT')
all_arg_human_df_ATC<-merge(arg_ACT,all_arg_human_df,'arg')
head(all_arg_human_df_ATC)
dim(all_arg_human_df_ATC)
all_arg_human_df_ATC_freq<-(as.data.frame(table(all_arg_human_df_ATC$ACT))) %>% 
  mutate(Perc=100*Freq/sum(Freq)) %>%
  mutate(group='human') %>%
  subset(select=c(Var1,Perc,group))
all_arg_human_df_ATC_freq

all_arg_res_rpkm_swine<-subset(all_arg_res_rpkm_with_region,Host == "Swine")
all_arg_res_rpkm_swine<-all_arg_res_rpkm_swine[,-c(2,3,4)]
rownames(all_arg_res_rpkm_swine)<-all_arg_res_rpkm_swine[,1]
all_arg_res_rpkm_swine<-all_arg_res_rpkm_swine[,-1]
dim(all_arg_res_rpkm_swine)
all_arg_res_rpkm_swine_clean<-all_arg_res_rpkm_swine[, colSums(all_arg_res_rpkm_swine != 0) > 0]
dim(all_arg_res_rpkm_swine_clean)
all_arg_swine_df<-data.frame(colnames(all_arg_res_rpkm_swine_clean))
colnames(all_arg_swine_df)<-'arg'
all_arg_swine_df_ATC<-merge(arg_ACT,all_arg_swine_df,'arg')
all_arg_swine_df_ATC_freq<-(as.data.frame(table(all_arg_swine_df_ATC$ACT))) %>% 
  mutate(Perc=100*Freq/sum(Freq))%>%
  mutate(group='swine') %>%
  subset(select=c(Var1,Perc,group))
all_arg_swine_df_ATC_freq

all_arg_res_rpkm_chicken<-subset(all_arg_res_rpkm_with_region,Host == "Chicken")
all_arg_res_rpkm_chicken<-all_arg_res_rpkm_chicken[,-c(2,3,4)]
rownames(all_arg_res_rpkm_chicken)<-all_arg_res_rpkm_chicken[,1]
all_arg_res_rpkm_chicken<-all_arg_res_rpkm_chicken[,-1]
dim(all_arg_res_rpkm_chicken)
all_arg_res_rpkm_chicken_clean<-all_arg_res_rpkm_chicken[, colSums(all_arg_res_rpkm_chicken != 0) > 0]
dim(all_arg_res_rpkm_chicken_clean)
all_arg_chicken_df<-data.frame(colnames(all_arg_res_rpkm_chicken_clean))
colnames(all_arg_chicken_df)<-'arg'
all_arg_chicken_df_ATC<-merge(arg_ACT,all_arg_chicken_df,'arg')
all_arg_chicken_df_ATC_freq<-(as.data.frame(table(all_arg_chicken_df_ATC$ACT))) %>% 
  mutate(Perc=100*Freq/sum(Freq))%>%
  mutate(group='chicken') %>%
  subset(select=c(Var1,Perc,group))
all_arg_chicken_df_ATC_freq

all_arg_df_ATC_freq<-rbind(all_arg_swine_df_ATC_freq,all_arg_human_df_ATC_freq,all_arg_chicken_df_ATC_freq)
all_arg_df_ATC_freq

library(viridis)
library(hrbrthemes)
library(scales)
colors=pal_npg('nrc')(10) ###library(scales)
show_col(colors)
pal_simpsons("springfield")(16)
show_col(pal_simpsons("springfield")(16))

fig3_p1<-ggplot(all_arg_df_ATC_freq, aes(fill=group, y=Perc, x=reorder(Var1,Perc))) + 
  geom_bar(position="stack",stat="identity",alpha=0.8) +
  #scale_fill_manual(values=c("#E64B35FF","#3C5488FF"))+
  scale_fill_manual(values=c(pal_npg(c("nrc"))(10)[c(1,4,3)]))+
  #scale_color_manual(values=c(pal_npg(c("nrc"))(10)[c(1,4)]))+
  #ggtitle("Classes of ARGs in both types of hosts") +
  theme_bw() + #theme_ipsum()+theme_minimal
  theme(axis.text.x = element_text(face='bold',family="Times",size = 22,angle = 90,vjust=0.6),
        axis.text.y=element_text(face='bold',family="Times",size = 22),
        axis.title = element_text(family="Times",face="bold",size = 24),
        legend.title = element_text(colour="black", family="Times",size=24, face="bold"),
        legend.text = element_text(colour="black", family="Times",size = 22))+ #axis.ticks.x = element_blank()
  xlab("")+
  scale_y_continuous(name = "Percentage of ARGs in each class(%)")+
  labs(fill = "Host")
fig3_p1

library(VennDiagram)
venn.plot<-venn.diagram(list(human=colnames(all_arg_res_rpkm_human_clean),
                             swine=colnames(all_arg_res_rpkm_swine_clean),
                             chicken=colnames(all_arg_res_rpkm_chicken_clean)),
                        filename = NULL,
                        #'China_HK_swine_venn.tiff',compression="lzw",
                        cex = 3,
                        cat.default.pos = "outer",
                        cat.cex = 3,
                        cat.dist = c(0.000, 0.000, 0.000),
                        cat.pos = c(-45, 45, 135),
                        fill = c(alpha("#E64B35FF",0.6),alpha("#3C5488FF",0.6),alpha("#00A087FF",0.6)),
                        #cat.default.pos = "text",
                        resolution = 300,units = 'px',lwd = 4) ### cat.default.pos = "outer" to let the caption outside

cowplot::plot_grid(venn.plot,p1,labels = c('a.', 'b.'),rel_widths=c(6,10))
grid.draw(venn.plot)
cowplot::plot_grid(venn.plot)

arg_tax_is_plasmids_1174<-read.table('./data/arg_tax_is_plasmids_1487.txt',
                                     stringsAsFactors = FALSE, sep='\t', head=F, quote = "#")
##arg_tax_is_plasmids_1174_IS<-subset(arg_tax_is_plasmids_1174,!is.na(arg_tax_is_plasmids_1174$V10))
dim(arg_tax_is_plasmids_1174)
sample_region
arg_tax_is_plasmids_1174_fil<-merge(sample_region,arg_tax_is_plasmids_1174,by.x='sample',by.y='V1')
arg_tax_is_plasmids_1174_fil<-merge(arg_ACT,arg_tax_is_plasmids_1174_fil,by.x='arg',by.y='V7')
head(arg_tax_is_plasmids_1174_fil)
dim(arg_tax_is_plasmids_1174_fil)
arg_tax_is_plasmids_1174_clean<-subset(arg_tax_is_plasmids_1174_fil,V9!='')
dim(arg_tax_is_plasmids_1174_clean)
arg_tax_is_plasmids_1174_unannotated<-subset(arg_tax_is_plasmids_1174_fil,is.na(V9))
dim(arg_tax_is_plasmids_1174_unannotated)

arg_tax_is_plasmids_1174_plasmids<-arg_tax_is_plasmids_1174_clean[str_detect(arg_tax_is_plasmids_1174_clean$V9,'plasmid'),] ## or arg_tax_is_plasmids_1174_plasmids[grep("plasmid", arg_tax_is_plasmids_1174_plasmids$V9), ]
head(arg_tax_is_plasmids_1174_plasmids)
dim(arg_tax_is_plasmids_1174_plasmids)
arg_tax_is_plasmids_1174_plasmids_df<-as.data.frame(dcast(arg_tax_is_plasmids_1174_plasmids[c('ACT','Host')],Host~ACT))
arg_tax_is_plasmids_1174_plasmids_df<-data.frame(arg_tax_is_plasmids_1174_plasmids_df,check.names = FALSE,row.names = arg_tax_is_plasmids_1174_plasmids_df$Host,stringsAsFactors = FALSE)
arg_tax_is_plasmids_1174_plasmids_df<-arg_tax_is_plasmids_1174_plasmids_df[,-1]
arg_tax_is_plasmids_1174_plasmids_df<-data.matrix(arg_tax_is_plasmids_1174_plasmids_df)
arg_tax_is_plasmids_1174_plasmids_df<-t(arg_tax_is_plasmids_1174_plasmids_df)
head(arg_tax_is_plasmids_1174_plasmids_df)
arg_tax_is_plasmids_1174_plasmids_df<-data.frame(arg_tax_is_plasmids_1174_plasmids_df)%>%
  rownames_to_column('ACT')%>%
  mutate(data.frame(arg_tax_is_plasmids_1174_plasmids_df),Human_plasmids=100*Human/sum(Human)) %>%
  mutate(Swine_plasmids=100*Swine/sum(Swine)) %>%
  mutate(Chicken_plasmids=100*Chicken/sum(Chicken))%>%
  column_to_rownames('ACT')
arg_tax_is_plasmids_1174_plasmids_df<-arg_tax_is_plasmids_1174_plasmids_df[,c('Human_plasmids','Swine_plasmids','Chicken_plasmids')]
sum(arg_tax_is_plasmids_1174_plasmids_df$Human)

#arg_tax_is_plasmids_1174_chromosome<-arg_tax_is_plasmids_1174_clean[str_detect(arg_tax_is_plasmids_1174_clean$V9,'chromosome'),]
arg_tax_is_plasmids_1174_chromosome<-arg_tax_is_plasmids_1174_clean[!str_detect(arg_tax_is_plasmids_1174_clean$V9,'plasmid'),]
arg_tax_is_plasmids_1174_chromosome<-arg_tax_is_plasmids_1174_chromosome[!str_detect(arg_tax_is_plasmids_1174_chromosome$V9,'unclassified.unclassified'),]
dim(arg_tax_is_plasmids_1174_chromosome)
head(arg_tax_is_plasmids_1174_chromosome)
arg_tax_is_plasmids_1174_chromosome_df<-as.data.frame(dcast(arg_tax_is_plasmids_1174_chromosome[c('ACT','Host')],Host~ACT))
arg_tax_is_plasmids_1174_chromosome_df<-data.frame(arg_tax_is_plasmids_1174_chromosome_df,check.names = FALSE,row.names = arg_tax_is_plasmids_1174_chromosome_df$Host,stringsAsFactors = FALSE)
arg_tax_is_plasmids_1174_chromosome_df<-arg_tax_is_plasmids_1174_chromosome_df[,-1]
arg_tax_is_plasmids_1174_chromosome_df<-data.matrix(arg_tax_is_plasmids_1174_chromosome_df)
arg_tax_is_plasmids_1174_chromosome_df<-t(arg_tax_is_plasmids_1174_chromosome_df)
head(arg_tax_is_plasmids_1174_chromosome_df)
arg_tax_is_plasmids_1174_chromosome_df<-data.frame(arg_tax_is_plasmids_1174_chromosome_df)%>%
  rownames_to_column('ACT')%>%
  mutate(data.frame(arg_tax_is_plasmids_1174_chromosome_df),Human_chromosome=100*Human/sum(Human)) %>%
  mutate(Swine_chromosome=100*Swine/sum(Swine))%>%
  mutate(Chicken_chromosome=100*Chicken/sum(Chicken))%>%
  column_to_rownames('ACT')
arg_tax_is_plasmids_1174_chromosome_df<-arg_tax_is_plasmids_1174_chromosome_df[,c('Human_chromosome','Swine_chromosome','Chicken_chromosome')]

arg_tax_is_plasmids_1174_unclassified<-arg_tax_is_plasmids_1174_clean[str_detect(arg_tax_is_plasmids_1174_clean$V9,c('unclassified.unclassified')),]
dim(arg_tax_is_plasmids_1174_unclassified)
head(arg_tax_is_plasmids_1174_unclassified)
arg_tax_is_plasmids_1174_unclassified_df<-as.data.frame(dcast(arg_tax_is_plasmids_1174_unclassified[c('ACT','Host')],Host~ACT))
arg_tax_is_plasmids_1174_unclassified_df<-data.frame(arg_tax_is_plasmids_1174_unclassified_df,check.names = FALSE,row.names = arg_tax_is_plasmids_1174_unclassified_df$Host,stringsAsFactors = FALSE)
arg_tax_is_plasmids_1174_unclassified_df<-arg_tax_is_plasmids_1174_unclassified_df[,-1]
arg_tax_is_plasmids_1174_unclassified_df<-data.matrix(arg_tax_is_plasmids_1174_unclassified_df)
arg_tax_is_plasmids_1174_unclassified_df<-t(arg_tax_is_plasmids_1174_unclassified_df)
head(arg_tax_is_plasmids_1174_unclassified_df)
dim(arg_tax_is_plasmids_1174_unclassified_df)
arg_tax_is_plasmids_1174_unclassified_df<-data.frame(arg_tax_is_plasmids_1174_unclassified_df)%>%
  rownames_to_column('ACT')%>%
  mutate(data.frame(arg_tax_is_plasmids_1174_unclassified_df),Human_unclassified=100*Human/sum(Human)) %>%
  mutate(Swine_unclassified=100*Swine/sum(Swine))%>%
  mutate(Chicken_unclassified=100*Chicken/sum(Chicken))%>%
  column_to_rownames('ACT')
arg_tax_is_plasmids_1174_unclassified_df<-arg_tax_is_plasmids_1174_unclassified_df[,c('Human_unclassified','Swine_unclassified','Chicken_unclassified')]
arg_tax_is_plasmids_1174_unclassified_df

arg_tax_is_MGE_1174_df<-merge(arg_tax_is_plasmids_1174_plasmids_df,arg_tax_is_plasmids_1174_chromosome_df,by='row.names')
arg_tax_is_MGE_1174_df<-merge(arg_tax_is_MGE_1174_df,arg_tax_is_plasmids_1174_unclassified_df,by.x='Row.names',by.y='row.names')
arg_tax_is_MGE_1174_df<-column_to_rownames(arg_tax_is_MGE_1174_df,'Row.names')
arg_tax_is_MGE_1174_df
library(pheatmap)
library(RColorBrewer)

fig3_p2<-pheatmap(as.matrix(arg_tax_is_MGE_1174_df),fontsize = 18,family="Times",
             fontfamily  = "Times",fontface = "bold",
             color=(colorRampPalette(brewer.pal(9, "Blues"))(1000)),
             display_numbers = FALSE,number_color='Red',
             cluster_cols=FALSE,cluster_rows = FALSE,
             angle_col = c("315"),
             silent=FALSE)
fig3_p2

cowplot::plot_grid(venn.plot,p1,p2$gtable,nrow = 3,ncol = 1,
                   #rel_widths=c(6,12),
                   rel_hight=c(2,2,2),
                   labels = c('a.', 'b.', 'c.', 'd.'),
                   label_size = 20)

library(circlize)
circos.clear()
arg_tax_is_plasmids_1174<-read.table('./data/arg_tax_is_plasmids_1487.txt', stringsAsFactors = FALSE, sep='\t', head=F, quote = "#")
##arg_tax_is_plasmids_1174_IS<-subset(arg_tax_is_plasmids_1174,!is.na(arg_tax_is_plasmids_1174$V10))
arg_tax_is_plasmids_1174_IS<-subset(arg_tax_is_plasmids_1174,V10!='')
head(arg_tax_is_plasmids_1174_IS)
dim(arg_tax_is_plasmids_1174_IS)

arg_sample<-arg_tax_is_plasmids_1174_IS[,c('V1','V7')]
head(arg_sample)

sample_region<-read.csv('./data/1487_sample.txt',header = F,sep='\t',row.names = 1,check.names = FALSE)
sample_region<-rownames_to_column(sample_region, var = "sample")
colnames(sample_region)<-c("sample",'region','type',"nation")
sample_region<-sample_region[,c(1,2)]
head(sample_region)
dim(sample_region)

arg_sample_c<-merge(sample_region,arg_sample,by.y ='V1',by.x = 'sample')
head(arg_sample_c)

arg_sample_c<-dcast(arg_sample_c,region~V7)
head(arg_sample_c)
dim(arg_sample_c)

arg_sample_c<-data.frame(arg_sample_c,check.names = FALSE,row.names = arg_sample_c$region,stringsAsFactors = FALSE)
arg_sample_c<-arg_sample_c[,-1]
arg_sample_c<-data.matrix(arg_sample_c)
arg_sample_t<-t(arg_sample_c)

arg_ACT<-read.csv('./data/arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
library(tibble)
arg_ACT<-rownames_to_column(arg_ACT, var = "arg")
colnames(arg_ACT)<-c("arg",'ACT')
arg_sample_t_class<-merge(arg_ACT,arg_sample_t,by.x = 'arg',by.y = 'row.names')
head(arg_sample_t_class)
dim(arg_sample_t_class)
unique(arg_sample_t_class$ACT)
library(dplyr)
library(tidyverse)
arg_sample_atc_sum<-apply(arg_sample_t_class[,3:38],2,function(x) tapply(x,arg_sample_t_class$ACT,sum))
arg_sample_atc_sum[is.na(arg_sample_atc_sum)] <- 0
colnames(arg_sample_atc_sum)
rownames(arg_sample_atc_sum)

colnames(arg_sample_atc_sum)<-c('AH','BeC','BeS','BuC','BuS','CaH','ChC','ChH','ChS','DC','DH','DS','FC','FH','FS','GC','GH',
                                'GS','HH','IcH','InH','ItC','ItH','ItS','JH','NC','NS','PH','PC','PS','SaH','SpC','SpH','SpS','SwH','UH')
rownames(arg_sample_atc_sum)<-c('Ami','Bla','Col','Fos','Mac','Nit','Phe','Qui','Rif','Sul','Tet','Tri','Van')

mycolors<-c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", 
            "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", 
            "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455")

grid.col = c(AH = "#771155", Bes = "#AA4488",Bus= "#CC99BB",CaH = "#114477",ChH = "#4477AA", ChS = "#4477AA",DH = "#77AADD",
             DS = "#77AADD", FH = "#117777", FS = "#117777", GH = "#44AAAA", GS = "#44AAAA", HH = "#77CCCC",
             IcH = "#117744", InH = "#44AA77", ItH = "#88CCAA", ItS = "#88CCAA", JH = "#777711", NS = "#AAAA44",
             PH = "#DDDD77", PS = "#DDDD77", SaH = "#774411", SpH = "#AA7744", SpS = "#AA7744", SwH = "#DDAA77",
             UH = "#771122",
             Ami = )

fig3d_p3<-chordDiagram(arg_sample_atc_sum,circos.axis(h='top'))

write.table(arg_sample_c,file = "arg_region_circle.txt",row.names = TRUE, col.names = TRUE,quote=FALSE,sep = '\t')

circos.par(gap.after = c(rep(2, nrow(arg_sample_atc_sum)-1), 10, rep(2, ncol(arg_sample_atc_sum)-1), 10))
chordDiagram(arg_sample_atc_sum, annotationTrack = "grid", transparency = 0.5,
             preAllocateTracks = list(track.height = 0.01),grid.col = grid.col)
#for(si in get.all.sector.index()) {
#  circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2)
#}
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  #circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3)
  #for(p in seq(0, 1, by = 0.25)) {
  #  circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim), p, facing = "clockwise",cex = 0.4, adj = c(0, 0.5), niceFacing = TRUE)
  # }#adj = c(0.5, -0.2)
  circos.text(mean(xlim), 1.4, cex=1.5,sector.name,niceFacing = TRUE)#facing = "clockwise", adj = c(0, 0.5),
}, bg.border = NA)

circos.clear()

p3

### network
suppressMessages(library(igraph))
arg_tax_is_plasmids_1174<-read.table('./data/arg_tax_is_plasmids_1487.txt',
                                     stringsAsFactors = FALSE, sep='\t', head=F, quote = "#")
##arg_tax_is_plasmids_1174_IS<-subset(arg_tax_is_plasmids_1174,!is.na(arg_tax_is_plasmids_1174$V10))
dim(arg_tax_is_plasmids_1174)
sample_region<-read.csv('./data/1487_sample.txt',header = F,sep='\t',row.names = 1,check.names = FALSE)
sample_region<-rownames_to_column(sample_region, var = "sample")
colnames(sample_region)<-c("sample",'region','type',"nation")
sample_region<-sample_region[,c(1,2)]
head(sample_region)
dim(ample_region)
arg_tax_is_plasmids_1174_fil<-merge(sample_region,arg_tax_is_plasmids_1174,by.x='sample',by.y='V1')
arg_tax_is_plasmids_1174_fil<-merge(arg_ACT,arg_tax_is_plasmids_1174_fil,by.x='arg',by.y='V7')
dim(sample_region)
dim(arg_ACT)
head(arg_tax_is_plasmids_1174_fil)
dim(arg_tax_is_plasmids_1174_fil)
arg_tax_is_plasmids_1174_IS<-subset(arg_tax_is_plasmids_1174_fil,V10!='')
head(arg_tax_is_plasmids_1174_IS)
dim(arg_tax_is_plasmids_1174_IS)
###### type is host: human, chicken, swine
arg_is_1174_clean<-select(arg_tax_is_plasmids_1174_IS,c('arg', 'type',rpkm='V6',IS='V10')) %>%
  unite(arg_type, c(arg, type), sep=':') #,rpkm
head(arg_is_1174_clean)
dim(arg_is_1174_clean)

#s <- strsplit(arg_is_1174_clean$IS, split = ",") #%>%unlist() %>%unnest()
s<-lapply(arg_is_1174_clean$IS, function (x) unique(unlist(strsplit(as.character(x), ","))))
s
length(s)
arg_is_1174_fil<-data.frame(arg_type = rep(arg_is_1174_clean$arg_type, sapply(s, length)), IS = unlist(s))
dim(arg_is_1174_fil)
head(arg_is_1174_fil)

ar.count <- arg_is_1174_fil%>%
  count(arg_type,IS)
data.frame(ar.count)
dim(ar.count)
head(ar.count)

cluster.count<-arg_is_1174_fil%>%
  group_by(arg_type) %>%   ## group by annotation
  mutate(cluster.size=length(arg_type)) %>% ungroup()
data.frame(cluster.count)

###### type is host: human, chicken, swine
arg_is_1174_clean<-select(arg_tax_is_plasmids_1174_IS,c('arg', 'type',rpkm='V6',IS='V10')) %>%
  unite(arg_type, c(arg, type,rpkm), sep=':') #
head(arg_is_1174_clean)
dim(arg_is_1174_clean)

#s <- strsplit(arg_is_1174_clean$IS, split = ",") #%>%unlist() %>%unnest()
s<-lapply(arg_is_1174_clean$IS, function (x) unique(unlist(strsplit(as.character(x), ","))))
s
length(s)
arg_is_1174_fil<-data.frame(arg_type = rep(arg_is_1174_clean$arg_type, sapply(s, length)), IS = unlist(s))
dim(arg_is_1174_fil)
head(arg_is_1174_fil)

arg_rpkm<-arg_is_1174_fil%>%
  mutate(rpkm=str_split_fixed(arg_type, ":", 3)[,3])%>%
  mutate(arg_type=str_remove(arg_type, ':[0-9].*+$'))%>%
  group_by(arg_type, IS) %>%
  summarise(arg_size=sum(as.numeric(rpkm)))
head(data.frame(arg_rpkm))
dim(arg_rpkm)

IS_rpkm<-arg_is_1174_fil%>%
  mutate(rpkm=str_split_fixed(arg_type, ":", 3)[,3])%>%
  mutate(arg_type=str_remove(arg_type, ':[0-9].*+$'))%>%
  group_by(IS) %>%
  summarise(IS_size=sum(as.numeric(rpkm)))
head(data.frame(IS_rpkm))
dim(IS_rpkm)

arg_is_1174_dat<-merge(ar.count,cluster.count,by=c(1,2))%>%
  merge(arg_rpkm,by=c(1,2))%>%
  filter(cluster.size>=5 ) %>%
  mutate(score=n/cluster.size) %>%
  filter(score>= 0.6 )
dim(arg_is_1174_dat)
head(arg_is_1174_dat)
arg_is_1174_dist<-distinct(arg_is_1174_dat)
dim(arg_is_1174_dist)
head(arg_is_1174_dist)
write.table(arg_is_1174_dat,'./arg_is_1174_dist.txt', quote=F, sep='\t', row.names = F, col.names = T)
write.table(arg_is_1174_dist,'./arg_is_1174_cytoscape_dist.tsv', quote=F, sep='\t', row.names = F, col.names = T)

arg_ACT<-read.csv('./data/arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
library(tibble)
arg_ACT<-rownames_to_column(arg_ACT, var = "arg")
colnames(arg_ACT)<-c("arg",'ACT')

arg_is_1174_dist_ACT<-arg_is_1174_dist%>%
  mutate(arg=str_split_fixed(arg_type, ":", 2)[,1])%>%
  merge(arg_ACT,'arg')%>%
  select(arg_type,IS,n,cluster.size,score,ACT,arg_size)%>%
  unite(arg_type,arg_type,ACT,arg_size,sep=':')
head(arg_is_1174_dist_ACT)

g <- graph_from_data_frame(arg_is_1174_dist_ACT, directed=FALSE)
V(g)$type <- 'NA'
V(g)$type[!(str_detect(V(g)$name, ':'))] <- 'IS'
V(g)$type[str_detect(V(g)$name, 'Human')] <- 'Human'
V(g)$type[str_detect(V(g)$name, 'Swine')] <- 'Swine'
#V(g)$type[str_detect(V(g)$name, 'Swine')] <- 'FoodAnimal'
V(g)$type[str_detect(V(g)$name, 'Chicken')] <- 'Chicken'
#V(g)$type[str_detect(V(g)$name, 'Chicken')] <- 'FoodAnimal'

V(g)$lab <- V(g)$name
V(g)$color <- ifelse(V(g)$type=='IS', 'IS', ((str_split_fixed((V(g)$name),":", 4)[,3])))
V(g)$color
V(g)$lab <- ifelse((str_detect(V(g)$name, ':')),((str_split_fixed((V(g)$name),":", 2)[,1])),(V(g)$name))
#V(g)$type <- str_detect(V(g)$name, 'Human')
#V(g)$color <- ifelse(V(g)$type, 'ARGs', 'IS')
V(g)$lab[str_detect(V(g)$lab, 'ISCco2')] <- 'Tn5541/Tn5543'
V(g)$lab
E(g)$weight <- arg_is_1174_dist_ACT$score
V(g)$size <- ifelse(V(g)$type=='IS', 10000, ((str_split_fixed((V(g)$name),":", 4)[,4])))
#c_scale <- colorRamp(c('red','yellow','cyan','blue'))
#E(g)$color<- apply(c_scale(E(g)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
#V(g)$size[V(g)$type] <- 20
#V(g)$lab[V(g)$type] <- paste0(tmp[V(g)$type, 1], "\n", "(",tmp[V(g)$type,2], ")")
V(g)$size

suppressMessages(library(igraph))
library(ggraph)
library(foreach)
library(stringr)
library(tidyr)
library(dplyr)
library(ggsci)
library(RColorBrewer)
library(scales)
library(tibble)
library(cowplot)
library(dplyr)
library(tidyverse)
library(ggplot2)
show_col(pal_npg(c("nrc"))(10))

### very good source: http://mr.schochastics.net/netVizR.html
## to use gradient color: https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/
g
ggraph(g, 'nicely') +  #'nicely'
  geom_edge_arc(aes(col=weight),#alpha=0.8,###aes(width=weight),col=weight>=edge_weight_thre, lty=weight>=edge_weight_thre
                curvature = 0.0)+#,
  #end_cap=circle(3, 'mm'), start_cap=circle(3, 'mm')) +
  geom_node_point(aes(shape=type,color=color,size=log2(as.numeric(size)))) + #, size=50, fill=color,
  geom_node_text(aes(label = lab), size=5, fontface=4,repel = TRUE,family="Times") + #'italic'
  #scale_edge_colour_brewer(palette = "Set1")+
  #scale_edge_color_gradient(low=pal_npg(c("nrc"))(10)[c(5)],high=pal_simpsons("springfield")(16)[9])+
  scale_edge_color_gradient(low=pal_simpsons("springfield")(16)[6],
                            high=pal_simpsons("springfield")(16)[9],name='Score')+ #(low="#00AFBB", high="#E7B800")+ #(low = 'yellow', high = 'cyan')
  #scale_edge_fill_manual(values=c(colorRampPalette(brewer.pal(9, "Reds"))(1000)))+
  #scale_fill_manual(values=c(colorRampPalette(brewer.pal(9,'Blues'))(25)))+
  #scale_color_gradient(low="blue", high="red")+
  #scale_fill_gradientn(colours = terrain.colors(10))+
  #scale_edge_color_manual(values=c(E(g)$color)) + #values=c('black','red'),col=colorRampPalette(brewer.pal(9, "Blues"))(1000)
  #scale_edge_width_continuous(range=c(1, 4)) + 
  #scale_edge_linetype_manual(values=c('dotdash','solid')) + 
  #scale_radius(range=c(1,5)) +
  scale_size(name='Relative abundance\nof ARGs(log2rpkm)')+
  scale_color_manual(values=c(pal_simpsons("springfield")(16)[-3]),name = 'Antibiotic Class') + #values=c(pal_npg(c("nrc"))(10)[c(1,2,4)])
  
  theme_void() + 
  scale_shape_manual(values=c(15,16,17,18),name='Host')+
  theme(legend.text=element_text(size=12,family="Times"),
        legend.title = element_text(size=14,family="Times"))

### to extract specific ARGs
head(arg_tax_is_plasmids_1174_fil)
dim(arg_tax_is_plasmids_1174_fil)
arg_tax_is_plasmids_1174_ant6<-subset(arg_tax_is_plasmids_1174_fil,arg=="ant(6)-Ia")
head(arg_tax_is_plasmids_1174_ant6)
dim(arg_tax_is_plasmids_1174_ant6)
write.table(arg_tax_is_plasmids_1174_ant6,'./arg_tax_is_plasmids_1174_ant6.txt',quote=F, sep='\t', row.names = F, col.names = T)

##### my code to plot trees with IS and ARGs
library(ggtree)
library(ggnewscale)
library(phytools)
library(ggsci)
library(scales)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(phyloseq)
library(stringr)

colors <- pal_npg('nrc')(10)
show_col(colors)

tr<-midpoint.root(read.tree('./data/arg_1487_ant6_800_cdhit_mafft.fasta.contree'))
ant6_meta<-read.table('./data/arg_1487_ant6_800_cdhit_host_count.txt',head=FALSE, sep='\t', stringsAsFactors = FALSE)
colnames(ant6_meta)<-c('laeles','Human','Swine','Chicken')
#ant6_meta_c<-ant6_meta%>%
#  mutate(log2Human=scale(Human))
#ant6_meta_c$log2Human[is.infinite(ant6_meta_c$log2Human)] <- 0
ant6_meta$laeles <- chartr("$", "_", ant6_meta$laeles) ### ant6_meta$laeles <- gsub("$", "_", ant6_meta$laeles)
ant6_meta<-data.frame(ant6_meta, row.names = 1)
ant6_meta$Swine[which(ant6_meta$Swine!=0)]<-1
ant6_meta$Human[which(ant6_meta$Human!=0)]<-1
ant6_meta$Chicken[which(ant6_meta$Chicken!=0)]<-1
ant6_meta$Human <- factor(ant6_meta$Human)
ant6_meta$Swine <- factor(ant6_meta$Swine)
ant6_meta$Chicken <- factor(ant6_meta$Chicken)
#ant6_meta <- mutate_if(ant6_meta,is.integer, as.factor)
dim(ant6_meta)
head(ant6_meta)
#rownames(ant6_meta)<-ant6_meta[,c(1)]
#ant6_meta<-ant6_meta[,-c(1)]

p<-ggtree(tr,lwd=1.5,open.angle = 10,layout="fan") 
p
p1_human<-gheatmap(p,ant6_meta['Human'],offset=0, width=.1,colnames = F,#offset=1.1
             colnames_angle=90, colnames_offset_y = .15,hjust=1.1)+ #,colnames_offset_x = .25,colnames_offset_y = .35
  scale_fill_manual(values=c('white',pal_npg(c("nrc"))(10)[1]))+
  #scale_fill_gradient(low = 'white',high = pal_npg(c("nrc"))(10)[1])+
  theme_tree(legend.title = element_blank())
p1_human
p2_tr <- p1_human + new_scale_fill()
p2_tr
p2_tr<-gheatmap(p2_tr,ant6_meta[c('Swine','Chicken')],offset=0.1, width=.2,colnames = F,#offset=1.1
             colnames_angle=90, colnames_offset_y = .35,hjust=1.1)+ #,colnames_offset_x = .25,colnames_offset_y = .35
  scale_fill_manual(values=c('white',pal_npg(c("nrc"))(10)[4]))+
  #scale_fill_gradient(low = 'white',high = pal_npg(c("nrc"))(10)[4])+
  theme_tree(legend.title = element_blank())
p2_tr


install.packages('gggenes')
library(gggenes)
library(ggplot2)
ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+
  theme_bw()


Tn4451_genes<-read.table('./Tn4451_genes_bed.txt',head=TRUE, sep='\t', stringsAsFactors = FALSE)

dummies <- make_alignment_dummies(Tn4451_genes,
                                  aes(xmin = start, xmax = end, y = molecule, id = gene,forward = direction),
                                  on = "tnpX")

ggplot(Tn4451_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene,forward = direction)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule,scales = "free",ncol = 1) + #scales = "free",
  #scale_fill_brewer(palette = "Set3")+
  scale_fill_simpsons()+
  geom_blank(data = dummies)+
  #scale_fill_npg()+
  theme_genes() ###theme_bw()

### to extract contigs of ARGs
head(arg_tax_is_plasmids_1174_fil)
dim(arg_tax_is_plasmids_1174_fil)
arg_tax_is_plasmids_1174_IS<-subset(arg_tax_is_plasmids_1174_fil,V10!='')
head(arg_tax_is_plasmids_1174_IS)
blaTEM_Tn2<-filter(arg_tax_is_plasmids_1174_IS, str_detect(arg_tax_is_plasmids_1174_IS$V10, 'Tn2'))%>%
  filter(str_detect(arg, 'blaTEM'))
head(blaTEM_Tn2)
dim(blaTEM_Tn2)

aph3III_ISCco2<-filter(arg_tax_is_plasmids_1174_IS, str_detect(arg_tax_is_plasmids_1174_IS$V10, 'ISCco2')) %>%
  filter(str_detect(arg, pattern=fixed("aph(3')-III"))) %>%
  select(sample,Region,contig=V2,tax=V8) %>%
  distinct() %>%
  mutate(length=str_split_fixed(contig, "_", 6)[,4],cov=str_split_fixed(contig, "_", 6)[,6])%>%
  filter(as.numeric(length)>=3000&as.numeric(cov)>=10)

aph3III_ISCco2<-filter(arg_tax_is_plasmids_1174_IS, str_detect(arg_tax_is_plasmids_1174_IS$V10, 'ISCco2')) %>%
  filter(str_detect(arg, pattern=fixed("aph(3')-III"))) %>%
  select(sample,Region,contig=V2,tax=V8) %>%
  distinct() %>%
  mutate(length=str_split_fixed(contig, "_", 6)[,4],cov=str_split_fixed(contig, "_", 6)[,6])%>%
  filter(str_detect(tax, pattern=fixed("Clostridioides")))

head(aph3III_ISCco2)
dim(aph3III_ISCco2)
hist(as.numeric(aph3III_ISCco2$length))

aac6_aph2_ISCco2<-filter(arg_tax_is_plasmids_1174_IS, str_detect(arg_tax_is_plasmids_1174_IS$V10, 'ISCco2'))%>%
  filter(str_detect(arg, pattern = fixed("aac(6')-aph(2'')")))

aac6Ia_ISCco2<-filter(arg_tax_is_plasmids_1174_IS, str_detect(arg_tax_is_plasmids_1174_IS$V10, 'ISCco2'))%>%
  filter(str_detect(arg, pattern = fixed("ant(6)-Ia")))

arg_ISCco2<-filter(arg_tax_is_plasmids_1174_IS, str_detect(arg_tax_is_plasmids_1174_IS$V10, 'ISCco2')) %>%
  #filter(str_detect(arg, pattern=fixed("aph(3')-III"))) %>%
  select(sample,Region,contig=V2,tax=V8) %>%
  distinct() %>%
  mutate(length=str_split_fixed(contig, "_", 6)[,4],cov=str_split_fixed(contig, "_", 6)[,6])#%>%
#filter(str_detect(tax, pattern=fixed("Clostridioides")))%>%
#filter((as.numeric(length)>=3000&as.numeric(cov)>=10)|length=='')

head(arg_ISCco2)
dim(arg_ISCco2)
write.table(arg_ISCco2,'./arg_1174_ISCco2.txt',quote=F, sep='\t', row.names = F, col.names = T)
unclassified_ISCco2<-filter(arg_ISCco2,tax=='root (taxid 1)' | tax == 'Bacteria (taxid 2)')
dim(unclassified_ISCco2)
write.table(unclassified_ISCco2,'./arg_1174_ISCco2_unclassified.txt',quote=F, sep='\t', row.names = F, col.names = T)

arg_ISCco2_tax<-data.frame(sort(table(arg_ISCco2$tax)))%>%
  mutate(Perc=100*Freq/sum(Freq))

others_perc<-sum(arg_ISCco2_tax$Perc[(arg_ISCco2_tax$Perc < 1.0)])
others_freq<-sum(arg_ISCco2_tax$Freq[(arg_ISCco2_tax$Perc < 1.0)])
arg_ISCco2_tax<-arg_ISCco2_tax%>%
  filter(Perc>1.0)%>%
  add_row(Var1 = 'others(<1%)',Freq = others_freq, Perc = others_perc)%>%
  mutate(species=str_extract(Var1, regex('[A-Z][a-z]+ [a-z]+')))
arg_ISCco2_tax$species[arg_ISCco2_tax$Var1 == 'root (taxid 1)']<-'Unclassified contigs'
arg_ISCco2_tax$species[arg_ISCco2_tax$Var1 == 'Bacteria (taxid 2)']<-'Unclassified Bacteria'
arg_ISCco2_tax$species[arg_ISCco2_tax$Var1 == 'others(<1%)']<-'others(<1%)'

ggplot(arg_ISCco2_tax, aes(x=species, y=Perc), lwd=2)+
  geom_bar(stat="identity", fill=pal_npg(c("nrc"))(10)[1],alpha=0.6)+
  labs(x='Tax',y='Percent of each tax')+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x=element_text(),
        axis.title.y=element_text(),
        axis.text.x=element_text(angle = 90),
        axis.text.y=element_text())

fig3g_p1<-ggplot(arg_ISCco2_tax, aes(x="",y=Perc,fill=reorder(species,-Perc)))+
  geom_bar(stat="identity",width=1)+#, fill=pal_npg(c("nrc"))(10)[1],alpha=0.6,
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values = pal_npg(c("nrc"))(10),name="Tn5541/Tn5543 Taxnomy")+
  theme(legend.text=element_text(size=16,family="Times",face = 4),
        legend.title = element_text(size=18,family="Times",face='bold'))
fig3g_p1
p1<-ggplot(arg_ISCco2_tax, aes(x="",y=Perc,fill=species))+
  geom_bar(stat="identity",width=1)+#, fill=pal_npg(c("nrc"))(10)[1],alpha=0.6,
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values = pal_npg(c("nrc"))(10),name="Tn5541/Tn5543 Taxnomy")+
  theme(legend.text=element_text(size=16,family="Times",face = 4),
        legend.title = element_text(size=18,family="Times",face='bold'))


arg_TnAs3<-filter(arg_tax_is_plasmids_1174_IS, str_detect(arg_tax_is_plasmids_1174_IS$V10, 'TnAs3')) %>%
  #filter(str_detect(arg, pattern=fixed("aph(3')-III"))) %>%
  select(sample,Region,contig=V2,tax=V8) %>%
  distinct() %>%
  mutate(length=str_split_fixed(contig, "_", 6)[,4],cov=str_split_fixed(contig, "_", 6)[,6]) #%>%
#filter(str_detect(tax, pattern=fixed("Aeromonas")))%>%
#filter(as.numeric(length)>=3000&as.numeric(cov)>=10)
head(arg_TnAs3)
dim(arg_TnAs3)
data.frame(sort(table(arg_TnAs3$tax)))
write.table(arg_TnAs3,'./arg_TnAs3_contigs.txt',quote=F, sep='\t', row.names = F, col.names = T)

arg_TnAs3_Aeromonas<-filter(arg_tax_is_plasmids_1174_IS, str_detect(arg_tax_is_plasmids_1174_IS$V10, 'TnAs3')) %>%
  #filter(str_detect(arg, pattern=fixed("aph(3')-III"))) %>%
  select(sample,Region,contig=V2,tax=V8) %>%
  distinct() %>%
  mutate(length=str_split_fixed(contig, "_", 6)[,4],cov=str_split_fixed(contig, "_", 6)[,6]) %>%
  filter(str_detect(tax, pattern=fixed("Aeromonas")))#%>%
#filter(as.numeric(length)>=3000&as.numeric(cov)>=10)
head(arg_TnAs3_Aeromonas)
dim(arg_TnAs3_Aeromonas)
write.table(arg_TnAs3_Aeromonas,'./arg_TnAs3_Aeromonas.txt',quote=F, sep='\t', row.names = F, col.names = T)

arg_TnAs3_tax<-data.frame(sort(table(arg_TnAs3$tax)))%>%
  mutate(Perc=100*Freq/sum(Freq))

others_perc<-sum(arg_TnAs3_tax$Perc[(arg_TnAs3_tax$Perc < 1.0)])
others_freq<-sum(arg_TnAs3_tax$Freq[(arg_TnAs3_tax$Perc < 1.0)])
arg_TnAs3_tax<-arg_TnAs3_tax%>%
  filter(Perc>1.0)%>%
  add_row(Var1 = 'others(<1%)',Freq = others_freq, Perc = others_perc)%>%
  mutate(species=str_extract(Var1, regex('[A-Z][a-z]+ [a-z]+')))
arg_TnAs3_tax$species[arg_TnAs3_tax$Var1 == 'root (taxid 1)']<-'Unclassified contigs'
arg_TnAs3_tax$species[arg_TnAs3_tax$Var1 == 'Bacteria (taxid 2)']<-'Unclassified Bacteria'
arg_TnAs3_tax$species[arg_TnAs3_tax$Var1 == 'others(<1%)']<-'others(<1%)'
arg_TnAs3_tax$species[arg_TnAs3_tax$Var1 == 'Enterobacterales (taxid 91347)']<-'Unclassified Enterobacterales'
arg_TnAs3_tax$species[arg_TnAs3_tax$Var1 == 'Gammaproteobacteria (taxid 1236)']<-'Unclassified Gammaproteobacteria'
arg_TnAs3_tax$species[arg_TnAs3_tax$Var1 == 'Enterobacteriaceae (taxid 543)']<-'Unclassified Enterobacteriaceae'
arg_TnAs3_tax$species[arg_TnAs3_tax$Var1 == 'Proteobacteria (taxid 1224)']<-'Unclassified Proteobacteria'

head(arg_TnAs3_tax)
fig3h_p2<-ggplot(arg_TnAs3_tax, aes(x="",y=Perc,fill=reorder(species,-Perc)))+
  geom_bar(stat="identity",width=1)+#, fill=pal_npg(c("nrc"))(10)[1],alpha=0.6,
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values = pal_simpsons("springfield")(16),name="TnAs3 Taxnomy")+
  theme(legend.text=element_text(size=16,family="Times",face = 4),
        legend.title = element_text(size=18,family="Times",face='bold'))
fig3h_p2
fig3h_p2<-ggplot(arg_TnAs3_tax, aes(x="",y=Perc,fill=species))+
  geom_bar(stat="identity",width=1)+#, fill=pal_npg(c("nrc"))(10)[1],alpha=0.6,
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_manual(values = pal_simpsons("springfield")(16),name="TnAs3 Taxnomy")+
  theme(legend.text=element_text(size=26,family="Times",face = 4),
        legend.title = element_text(size=28,family="Times",face='bold'))
