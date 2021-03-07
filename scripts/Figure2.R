##### code to all 1174 samples with 655 ARGs #####
library('gplots')
library(vegan)
library(data.table)
library(ggpubr)
library(readr)
library(dplyr)
library(readxl)
library('tibble') #### rownames_to_column
library(ggplot2)
library(reshape2)
library(ggsci)
library(ComplexHeatmap) ### to heatmap  
### https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap
### https://www.jianshu.com/p/86ae39a227f4
library(RColorBrewer) ### to brewer.pal
library(tidyverse)

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

sample_region<-read.csv('data/1487_sample.txt',header = F,sep='\t',row.names = 1,check.names = FALSE)
sample_region<-rownames_to_column(sample_region, var = "sample")
colnames(sample_region)<-c("sample",'Type','Host',"Region")
head(sample_region)
dim(sample_region)
#write.csv(sample_region,file = "sample_region.csv",row.names = TRUE, col.names = TRUE,quote=FALSE)

all_arg_res_rpkm_with_region<-merge(sample_region,all_arg_res_rpkm_t,'sample')
dim(all_arg_res_rpkm_with_region)
head(all_arg_res_rpkm_with_region)
all_arg_res_rpkm_with_region$sum<-rowSums(all_arg_res_rpkm_with_region[,5:867])
all_arg_res_rpkm_with_region$log2sum<-log2(all_arg_res_rpkm_with_region$sum)
levels(all_arg_res_rpkm_with_region$Region)

p1<-ggboxplot(all_arg_res_rpkm_with_region, x="Host", y="log2sum",
              add = "jitter", color = "Host",
              ggtheme = theme_bw(),
              palette = pal_npg(c("nrc"))(10)[c(1,4,3)])+
  stat_compare_means(method = "wilcox.test", label.y = 18)
p1<-ggpar(p1,legend="none",
          xlab = "",
          ylab = "Abundance of ARGs(log2(RPKM))",
          font.xtickslab = c(12, "bold"),
          font.y = c(14, "bold"),
          font.family = "Times Roman")
p1
p1<-ggplot(all_arg_res_rpkm_with_region, aes(x=Host, y=log2sum))+#fill=Host
  geom_boxplot(alpha=0.7,notch=TRUE) +
  #facet_wrap(~type, scale="free")+
  #facet_grid(. ~ type)+
  scale_y_continuous(name = "Abundance of ARGs(log2(RPKM))",
                     breaks = seq(0, 20, 5),
                     limits=c(0, 20)) +
  scale_x_discrete(name = "Host") +
  ggtitle("Relative abundance of ARGs in each type of hosts") +
  #theme(panel.grid.minor = element_blank(),panel.background = element_blank()) +
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 14),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        legend.position = "none") +
  #scale_fill_manual(values=c(colorRampPalette(brewer.pal(9,'Blues'))(25)))+
  scale_fill_manual(values=pal_npg(c("nrc"))(10)[c(1,4,3)])+
  #scale_color_gradient(low="blue", high="red")+
  #scale_fill_grey()+
  #scale_fill_brewer(palette = "Accent") +
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(fill = "Type of samples")
  
p1
install.packages("ggstatsplot")
install.packages("statsExpressions")
library(ggstatsplot)
library(statsExpressions)

my_pal = pal_npg(c("nrc"))(10)


p1<-ggstatsplot::ggbetweenstats(
  data = all_arg_res_rpkm_with_region,
  x=Host, 
  y=log2sum,
  notch = TRUE, # show notched box plot
  mean.plotting = TRUE, # whether mean for each group is to be displayed
  mean.ci = TRUE, # whether to display confidence interval for means
  mean.label.size = 6, # size of the label for mean
  type = "parametric", # which type of test is to be run
  k = 3, # number of decimal places for statistical results
  outlier.tagging = TRUE, # whether outliers need to be tagged
  outlier.label = Type, # variable to be used for the outlier tag
  outlier.label.color = "darkgreen", # changing the color for the text label
  #outlier.label.size = 4,
  xlab = "", # label for the x-axis variable
  ylab = "Abundance of ARGs(log2(RPKM))", # label for the y-axis variable
  title = "", # title text for the plot
  ggtheme = ggthemes::theme_few(), # choosing a different theme,theme_fivethirtyeight()
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  #package = "ggsci", # package from which color palette is to be taken,"wesanderson"
  #palette = "nrc_npg", # choosing a different color palette,"Darjeeling1"
  messages = FALSE
)+ggplot2::scale_color_manual(values = pal_npg(c("nrc"))(10)[c(1,4,3)])+
  ggplot2::theme(plot.title = element_text(size = 12, face = "bold",family = "Times"),
        text = element_text(size = 12,family = "Times"),
        axis.title = element_text(face="bold",family = "Times"),
        axis.text.x=element_text(size = 12,family = "Times"),
        legend.position = "none")
  

p1

p2<-ggplot(all_arg_res_rpkm_with_region[-(which(all_arg_res_rpkm_with_region$Type=='IcelandHuman')),], 
           aes(x = Type, y = log2sum, fill = Host)) + #fill = region
  geom_boxplot(alpha=0.7) +
  facet_wrap(~Host, scale="free")+
  #facet_grid(. ~ type)+
  scale_y_continuous(name = "Abundance of ARGs(log2(RPKM))",
                     breaks = seq(0, 20, 5),
                     limits=c(0, 20)) +
  #scale_x_discrete(name = "Type of samples") +
  #ggtitle("Relative abundance of ARGs in each type of samples") +
  #theme(panel.grid.minor = element_blank(),panel.background = element_blank()) +
  theme_bw()+
  #scale_fill_manual(values=c(colorRampPalette(brewer.pal(9,'Blues'))(25)))+
  scale_fill_manual(values=pal_npg(c("nrc"))(10)[c(1,4,3)])+
  #scale_color_gradient(low="blue", high="red")+
  #scale_fill_grey()+
  #scale_fill_brewer(palette = "Accent") +
  #labs(fill = "Type of samples")+
  theme(plot.title = element_text(size = 24, face = "bold",family = "Times"),
      text = element_text(size = 24,family = "Times"),
      axis.title = element_text(face="bold",family = "Times"),
      axis.text.x=element_text(size = 24,family = "Times",angle = 90,hjust = 1.0,face="bold"),
      axis.title.x=element_blank(),
      strip.text.x = element_text(size = 28,face = "bold", family = "Times"),
      legend.position = "none")
p2<-p2+ geom_jitter(shape=16, position=position_jitter(0.2))
p2


####heatmap ACT with regions
# https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/

library(pheatmap)
library(grid)
data.frame(colnames(all_arg_res_rpkm_with_region)[5:867])

arg_ACT<-read.csv('./data/arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
arg_ACT<-rownames_to_column(arg_ACT, var = "arg")
colnames(arg_ACT)<-c("arg",'ACT')
dim(all_arg_res_rpkm_1487)
all_arg_res_rpkm_ATC<-rownames_to_column(all_arg_res_rpkm_1174,var='arg')
all_arg_res_rpkm_ATC<-merge(arg_ACT,all_arg_res_rpkm_ATC,'arg')
head(all_arg_res_rpkm_ATC)
dim(all_arg_res_rpkm_ATC)
all_arg_res_rpkm_ATC_sum<-apply(all_arg_res_rpkm_ATC[,3:1176],2,function(x) tapply(x,all_arg_res_rpkm_ATC$ACT,sum))
dim(all_arg_res_rpkm_ATC_sum)
head(all_arg_res_rpkm_ATC_sum)
dim(all_arg_res_rpkm_ATC_sum[,colSums(all_arg_res_rpkm_ATC_sum[])>0])
dim(all_arg_res_rpkm_ATC_sum[rowSums(all_arg_res_rpkm_ATC_sum[])>0,])

all_arg_res_rpkm_ATC_sum_log2<-log2(all_arg_res_rpkm_ATC_sum[,])
is.infinite(all_arg_res_rpkm_ATC_sum_log2) %>% table()
is.na(all_arg_res_rpkm_ATC_sum_log2) %>% table()
is.nan(all_arg_res_rpkm_ATC_sum_log2) %>% table()

all_arg_res_rpkm_ATC_sum_log2[is.infinite(all_arg_res_rpkm_ATC_sum_log2)] <- 0
all_arg_res_rpkm_ATC_sum_log2[is.na(all_arg_res_rpkm_ATC_sum_log2)] <- 0

dim(all_arg_res_rpkm_ATC_sum)
all_arg_res_rpkm_ATC_sum_log2[is.infinite(all_arg_res_rpkm_ATC_sum_log2)]


atc_annot_region<-data.frame(region=sample_region$Region,sample_id=colnames(all_arg_res_rpkm_ATC_sum))
dim(atc_annot_region)
annotation_col = data.frame(Region=atc_annot_region$region,row.names = atc_annot_region$sample_id)
#newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotation_col$Region))))
#mycolors <- newCols(length(unique(annotation_col$Region)))
mycolors <- colorRampPalette(brewer.pal(9,'Blues'))(length(unique(annotation_col$Region)))
names(mycolors) <- unique(annotation_col$Region)
mycolors <- list(Region = mycolors)
mycolors
is.na(as.matrix(scale(all_arg_res_rpkm_ATC_sum_log2))) %>% table()

all_arg_res_rpkm_ATC_sum_log2_fil <- all_arg_res_rpkm_ATC_sum_log2[apply(all_arg_res_rpkm_ATC_sum_log2, 1, function(x) sd(x)!=0),]
dim(all_arg_res_rpkm_ATC_sum_log2_fil)
all_arg_res_rpkm_ATC_sum_fil <- all_arg_res_rpkm_ATC_sum[apply(all_arg_res_rpkm_ATC_sum, 1, function(x) sd(x)!=0),]
all_arg_res_rpkm_ATC_sum_fil <- all_arg_res_rpkm_ATC_sum[,apply(all_arg_res_rpkm_ATC_sum, 2, function(x) sd(x)!=0)]
dim(all_arg_res_rpkm_ATC_sum_fil)
is.na(as.matrix(standardize(all_arg_res_rpkm_ATC_sum))) %>% table()
is.na(as.matrix(scale(all_arg_res_rpkm_ATC[,3:1489]))) %>% table() 

test_m<-as.matrix(scale(all_arg_res_rpkm_ATC_sum_log2))
test_m[is.na(test_m)] <- 0
test_m
p3<-pheatmap(as.matrix(scale(all_arg_res_rpkm_ATC_sum_log2)),
             annotation_col = annotation_col,show_colnames=FALSE,
             fontsize=8,
             color = colorRampPalette(c("#4DBBD5FF","white", "#E64B35FF"))(80))
             #,
             #,
             #cellheight = 0.5,
             #cellwidth = 0.5,
             #color=colorRampPalette(brewer.pal(9,'Reds'))(60))
             
             #annotation_colors = mycolors,silent = TRUE)
p3

p3<-pheatmap(as.matrix(scale(all_arg_res_rpkm_ATC_sum_log2)),
         annotation_col = annotation_col,show_colnames=FALSE,
         fontsize=8,
         #cellheight = 0.5,
         #cellwidth = 0.5,
         #color=colorRampPalette(brewer.pal(9,'Reds'))(60))
         color = colorRampPalette(c("#4DBBD5FF","white", "#E64B35FF"))(80),
         annotation_colors = mycolors,silent = TRUE)
         #annotation_colors=colorRampPalette(brewer.pal(9,'Blues'))(26))
p3

library(randomcoloR)
annotation_cols <- sample_region[,-c(2)]
rownames(annotation_cols)<-annotation_cols[,1]
mycolors<-c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", 
  "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", 
  "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455")
#mycolors<-c(distinctColorPalette(20))
names(mycolors) <- unique(annotation_cols$Region)

annotation_colors<-list(Host=c(Human=pal_npg(c("nrc"))(10)[1],Swine=pal_npg(c("nrc"))(10)[3],
                        Chicken=pal_npg(c("nrc"))(10)[4]),
                        Region=mycolors)

p3<-pheatmap(scale(all_arg_res_rpkm_ATC_sum_log2),#scale = "none",
             annotation_col = annotation_cols[,-c(1)],show_colnames=FALSE,
             fontsize=24,
             fontfamily  = "Times",
             legend = FALSE,
             annotation_legend = FALSE,
             #cellheight = 0.5,
             #cellwidth = 0.5,
             #color=colorRampPalette(brewer.pal(9,'Reds'))(60))
             color = colorRampPalette(c("#4DBBD5FF","white", "#E64B35FF"))(80),
             annotation_colors = annotation_colors)

(p3$gtable)$grobs[[4]]$label

annotation_cols[,-c(1)]
### to generate p4
library(tibble)
library(ComplexHeatmap) ### to heatmap  
### https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap
library(RColorBrewer) ### to brewer.pal
library(pheatmap)
library(gplots)
library(magrittr)
library(ggpubr)

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
#write.csv(sample_region,file = "sample_region.csv",row.names = TRUE, col.names = TRUE,quote=FALSE)

all_arg_res_rpkm_with_region<-merge(sample_region,all_arg_res_rpkm_t,'sample')
dim(all_arg_res_rpkm_with_region)
head(all_arg_res_rpkm_with_region)
all_arg_res_rpkm_with_region[is.na(all_arg_res_rpkm_with_region)] <- 0
all_arg_res_rpkm_with_region$sum<-rowSums(all_arg_res_rpkm_with_region[,5:867])
all_arg_res_rpkm_with_region$log2sum<-log2(all_arg_res_rpkm_with_region$sum)
dim(all_arg_res_rpkm_with_region)
levels(all_arg_res_rpkm_with_region$region)

### to draw heatmap for stastics significant difference ARGs between human and swine
rowMedians <- function(x) apply(x, 1, median) ## aux function
all_arg_res_rpkm_with_region_clean<-select(all_arg_res_rpkm_with_region,-c(sum,log2sum))
all_arg_res_rpkm_with_region.ordered<-all_arg_res_rpkm_with_region_clean[order(all_arg_res_rpkm_with_region_clean$type),]
all_arg_res_rpkm.ordered<-all_arg_res_rpkm_with_region.ordered[,-c(2,3,4)]
rownames(all_arg_res_rpkm.ordered)<-all_arg_res_rpkm.ordered[,1]
all_arg_res_rpkm.ordered<-all_arg_res_rpkm.ordered[,-1]
#all_arg_res_rpkm.ordered<-data.frame(all_arg_res_rpkm.ordered,check.names = FALSE)
all_arg_res_rpkm.ordered<-t(all_arg_res_rpkm.ordered)
medianChicken <- rowMedians(all_arg_res_rpkm.ordered[,1:313])
medianHuman <- rowMedians(all_arg_res_rpkm.ordered[,314:1010])
medianSwine <- rowMedians(all_arg_res_rpkm.ordered[,1011:1487])
data1 <- all_arg_res_rpkm.ordered[which(medianHuman > 0.5 | medianSwine > 0.5 | medianChicken > 0.5),]
g <- factor(rep(1:3, c(313, 697, 477)), labels = c("FA","H","FA"))
dim(data1)
rownames(data1)
p.values <- sapply(seq(dim(data1)[1]), function(x) wilcox.test(as.numeric(data1[x,(g=="H")]),
                                                               as.numeric(data1[x,g=="FA"]))$p.value)
q.values <- p.adjust(p.values, method="fdr")
p.values
q.values
g
data2<-data1
colnames(data2)<-g
data2t<-t(data2)%>%data.frame(check.names = T)
data2t<-rownames_to_column(data2t, var = "type")
kruskal.test(weight ~ group, data = data1)
pvalueD=p.adjust(sapply(colnames(data2t[,-1]), function(x) kruskal.test(as.formula(paste0(x, "~type")),data=data2t)$p.value),method="fdr")
pvalueD

#filter data to plot heatmap
data2 <- data1[which(q.values<0.01),]
dim(data2)
data2Human <- data2[rowMedians(data2[,1:697])>rowMedians(data2[,698:1174]),]
data2Human <- data2[(rowMedians(data2[,314:1010])>rowMedians(data2[,1:313]) & rowMedians(data2[,314:1010])>rowMedians(data2[,1011:1487])),]
data2Swine <- data2[rowMedians(data2[,1:697])<rowMedians(data2[,698:1174]),]
data2FA <- data2[(rowMedians(data2[,314:1010])<rowMedians(data2[,1:313]) & rowMedians(data2[,314:1010])<rowMedians(data2[,1011:1487])),]

data2VHuman <- cbind(rowMedians(data2Human[,1:697]),rowMedians(data2Human[,698:1174]))
data2VHuman <- cbind(rowMedians(data2Human[,314:1010]),rowMedians(cbind(data2Human[,1:313],data2Human[,1011:1487])))
data2VSwine <- cbind(rowMedians(data2Swine[,1:697]),rowMedians(data2Swine[,698:1174]))
data2VFA <- cbind(rowMedians(data2FA[,314:1010]),rowMedians(cbind(data2FA[,1:313],data2FA[,1011:1487])))

data2all <- rbind(data2VHuman,data2VSwine)
data2all <- rbind(data2VHuman,data2VFA)
rownames(data2all) <- rownames(rbind(data2Human,data2Swine))
rownames(data2all) <- rownames(rbind(data2Human,data2VFA))
colnames(data2all) <- c('Human', 'FoodAnimal')
head(data2all)
dim(data2all)

arg_ACT<-read.csv('data/arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
library(tibble)
arg_ACT<-rownames_to_column(arg_ACT, var = "arg")
colnames(arg_ACT)<-c("arg",'ACT')
arg_diff_all<-merge(arg_ACT,data2all,by.x = 'arg',by.y = 'row.names')
head(arg_diff_all)
dim(arg_diff_all)

arg_order<-data.frame((p3$gtable)$grobs[[4]]$label)
arg_order<-rownames_to_column(arg_order)
colnames(arg_order)<-c("order","ACT")
dim(arg_diff_all_ordered)
arg_diff_all<-merge(arg_order,arg_diff_all,by="ACT")
arg_diff_all_ordered<-arg_diff_all[order(as.numeric(arg_diff_all$order)),]
head(arg_diff_all_ordered)

row_anno_arg<-arg_diff_all_ordered%>%
  select(arg,ACT)
head(row_anno_arg)
#rownames(row_anno_arg)<-row_anno_arg[,1]
row_anno_arg<-subset(row_anno_arg, select = -c(arg) )
pal_npg(c("nrc"))(10)[1:8]
mycolors<-pal_npg(c("nrc"))(10)[c(1:6)]
names(mycolors) <- unique(row_anno_arg$ACT)

arg_diff_all_ordered<-arg_diff_all_ordered[,-c(1,2)]
rownames(arg_diff_all_ordered)<-arg_diff_all_ordered[,1]
arg_diff_all_ordered<-arg_diff_all_ordered[,-1]
head(arg_diff_all_ordered)
row_anno_arg<-row_anno_arg%>%
  select(ACT)
col_anno_host<-as.data.frame(group)

annotation_colors_p4<-list(Host=c(Human=pal_npg(c("nrc"))(10)[1],FoodAnimal=pal_npg(c("nrc"))(10)[4]),
                        ACT=mycolors)

p4<-pheatmap(log(arg_diff_all_ordered+1,2),
         show_colnames=T,#fontface="italic",
         color = colorRampPalette(brewer.pal(9,'Blues'))(80),
         cluster_rows = F, annotation_row = row_anno_arg,
         cluster_cols = F,
         cellwidth = 15,
         #cellheight = 15,
         fontsize = 12, fontsize_row = 10, fontsize_col = 18,
         fontfamily = 'Times',annotation_colors = annotation_colors_p4,
         border_color="black",silent=FALSE)
p4

dist.mat <- vegdist(t(all_arg_res_rpkm_1487),na.rm = T)
cmds <- cmdscale(dist.mat, k=3, eig=TRUE)
eigen <- cmds$eig / sum(cmds$eig) * 100
dat.merged <- (merge(cmds$points, sample_region, by.x=0, by.y="sample", all.x=TRUE))
head(dat.merged)
dat.merged$region <- relevel(factor(dat.merged$region), 'HKHuman')
## PCoA
p5<-ggplot(dat.merged, aes(x=V1, y=V2, col=Region,shape=Host), lwd=2) +
  geom_density_2d(aes(x=V1, y=V2), inherit.aes = FALSE, col='grey', lwd=1) + 
  geom_point(size=4.5) +
  labs(x=paste0('MDS2 (',round(eigen[1], 1),'%)'),
       y=paste0('MDS1 (',round(eigen[2], 1),'%)')) +
  scale_color_manual(name="Region",values = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", 
                                "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", 
                                "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455"))+
  #scale_color_manual(values=pal_npg(c("nrc"))(10)[c(1,5,7,2,3,10,4,6,8,9,1,5,7,2,3,10,4,6,8,9,1,5,7,2,3,10)]) + 
  scale_shape_discrete(name = "Host")+
#  scale_color_discrete(name = "Region")+
  theme(plot.title = element_text(size = 14, face = "bold",family = "Times"),
        text = element_text(size = 14,family = "Times"),
        axis.title = element_text(face="bold",family = "Times"),
        axis.text.x=element_text(size = 14,family = "Times",angle = 45,hjust = 1.0),
        #axis.title.x=element_blank(),
        legend.position = "bottom")
#  theme(legend.title = element_blank())
p5

##read data
all_arg_res_rpkm <- read.csv("../../microbiome/data/resfinder_all/1174all_arg_res_rpkm.txt",sep="\t", 
                             header = T, row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
all_arg_res_rpkm_t<-t(all_arg_res_rpkm_1487)
all_arg_res_rpkm_t<-data.frame(all_arg_res_rpkm_t,check.names = FALSE)
all_arg_res_rpkm_t<-rownames_to_column(all_arg_res_rpkm_t, var = "sample") #### tibble is needed
sample_region<-read.csv('1174_sample.txt',header = F,sep='\t',row.names = 1,check.names = FALSE)
sample_region<-rownames_to_column(sample_region, var = "sample")
colnames(sample_region)<-c("sample",'region','type','nation')
all_arg_res_rpkm_with_region<-merge(sample_region,all_arg_res_rpkm_t,'sample')
levels(all_arg_res_rpkm_with_region$nation)
list.countries=c("China","Denmark","France","German","HK","Italy","Spain")
all_arg_res_rpkm_pairs<-all_arg_res_rpkm_with_region[which(all_arg_res_rpkm_with_region$nation %in% list.countries),]
metadat_select=data.frame(select(all_arg_res_rpkm_pairs, sample,type,nation),row.names = "sample")
all_arg_res_rpkm_pairs_fil<-all_arg_res_rpkm_pairs[,-c(2,3,4)]
rownames(all_arg_res_rpkm_pairs_fil)<-all_arg_res_rpkm_pairs_fil[,1]
all_arg_res_rpkm_pairs_fil<-all_arg_res_rpkm_pairs_fil[,-1]

dist.mat <- vegdist(all_arg_res_rpkm_pairs_fil,na.rm = T)
cmds <- cmdscale(dist.mat, eig=TRUE)
eigen <- cmds$eig / sum(cmds$eig) * 100
dat.merged <- (merge(cmds$points, metadat_select, by.x=0, by.y=0, all.x=TRUE))
dat.merged.m<-dat.merged%>%mutate(col=nation)
dat.merged.m$nation[dat.merged.m$nation == "HK"] <- "China"

p6<-ggplot(dat.merged.m, aes(x=V1, y=V2, col=col, shape=type), lwd=2) +
  #geom_curve(data=plot.dat, aes(x=V1.x,y=V2.x,xend=V1.y, yend=V2.y),
  #           arrow = arrow(length = unit(0.02, "npc")), lwd=1, alpha=0.5,
  #           inherit.aes = FALSE) +
  geom_point(size=3, alpha=0.9) + coord_cartesian(ylim = c(-0.5, 0.5), xlim = c(-0.75,0.5)) +
  labs(x=paste0('MDS1 (',round(eigen[1], 1),'%)'),
       y=paste0('MDS2 (',round(eigen[2], 1),'%)')) +
  scale_shape_manual(values=c(17,19,15)) + 
  scale_color_manual(values=pal_npg(c("nrc"))(10)[c(1,5,7,2,3,10,4)]) +
  facet_wrap(~nation,dir="v") + 
  guides(color=guide_legend(title='Region'), shape=guide_legend('Host'))+
  theme(plot.title = element_text(size = 14, face = "bold",family = "Times"),
        text = element_text(size = 14,family = "Times"),
        axis.title = element_text(face="bold",family = "Times"),
        axis.text.x=element_text(size = 14,family = "Times",angle = 45,hjust = 1.0),
        #axis.title.x=element_blank(),
        strip.text.x = element_text(size = 16),
        legend.position = "right")
dim(dat.merged)
p6

library(cowplot)
top_row <- plot_grid(p1, p2, labels = c('a.', 'b.'), label_size = 25,rel_widths=c(40,16))
left_col<-plot_grid(top_row, p3$gtable,labels = c('', 'c.'), label_size = 25, ncol = 1,
                    rel_heights = c(2,2))
left_col<-plot_grid(top_row, p3$gtable,labels = c('', 'c.'), label_size = 25, ncol = 1,
                    rel_heights = c(2,2))
upper_row <- plot_grid(left_col, p4$gtable,labels = c('','d.'), label_size = 25, nrow = 1,rel_widths=c(18,4))
lower_row <-plot_grid(p5,p6,nrow = 1,ncol = 2,labels = c('e.','f.'), label_size = 25,
                      rel_widths=c(12,8))
#left_col
plot_grid(upper_row, lower_row, nrow = 2,ncol=1,rel_heights = c(4,2))

ggsave("./microbiota_new_figures/Figs/Fig2.pdf", height = 12, width = 15)
