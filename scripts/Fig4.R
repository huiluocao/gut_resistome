
install.packages('VennDiagram')
library(VennDiagram)

ChinaHuman_all_arg_res_rpkm<-subset(all_arg_res_rpkm_with_region,region=='ChinaHuman')
ChinaHuman_all_arg_res_rpkm_clean<-ChinaHuman_all_arg_res_rpkm[, colSums(ChinaHuman_all_arg_res_rpkm != 0) > 0]
ChinaSwine_all_arg_res_rpkm<-subset(all_arg_res_rpkm_with_region,region=='ChinaSwine')
ChinaSwine_all_arg_res_rpkm_clean<-ChinaSwine_all_arg_res_rpkm[, colSums(ChinaSwine_all_arg_res_rpkm != 0) > 0]
HKHuman_all_arg_res_rpkm<-subset(all_arg_res_rpkm_with_region,region=='HKHuman')
HKHuman_all_arg_res_rpkm_clean<-HKHuman_all_arg_res_rpkm[, colSums(HKHuman_all_arg_res_rpkm != 0) > 0]

venn.plot<-venn.diagram(list(HKH=colnames(HKHuman_all_arg_res_rpkm_clean)[5:314],
                             CNH=colnames(ChinaHuman_all_arg_res_rpkm_clean)[5:224],
                             CNS=colnames(ChinaSwine_all_arg_res_rpkm_clean)[5:329]),
                        filename = NULL,
                        #'China_HK_swine_venn.tiff',compression="lzw",
                        cex = 3,cat.cex = 3,cat.dist = c(0.08, 0.08, 0.08),
                        fill = c(alpha("#E64B35FF",0.4),alpha("#4DBBD5FF",0.4),alpha("#3C5488FF",0.4)),
                        cat.default.pos = "text",resolution = 300,units = 'px',lwd = 2) ### cat.default.pos = "outer" to let the caption outside

cowplot::plot_grid(venn.plot)

cowplot::plot_grid(venn.plot,p1,labels = c('a.', 'b.'))
grid.draw(venn.plot)

all_arg_sum_rpkm_region<-all_arg_res_rpkm_with_region[-c(1,3:4)] 
all_arg_sum_rpkm_HK_CN<-all_arg_sum_rpkm_region %>%  #dt[, !c("x", "y")],df[, !(colnames(df) %in% c("x","bar","foo"))],within(df, rm(x, y)),df <- subset(df, select = -c(a, c))
  group_by(region) %>%
  summarise_each(funs(sum))%>%
  data.frame(check.names = F,stringsAsFactors = F)%>%
  subset(region %in% c('ChinaHuman','HKHuman','ChinaSwine'))
all_arg_sum_rpkm_HK_CN<-all_arg_sum_rpkm_HK_CN[, colSums(all_arg_sum_rpkm_HK_CN != 0) > 0]
head(all_arg_sum_rpkm_HK_CN)
rownames(all_arg_sum_rpkm_HK_CN) <- all_arg_sum_rpkm_HK_CN[,1]
all_arg_sum_rpkm_HK_CN.t<-t(all_arg_sum_rpkm_HK_CN[-1])
all_arg_sum_rpkm_HKH_CNH_CNS<-merge(arg_ACT,all_arg_sum_rpkm_HK_CN.t,by.x='arg',by.y=0)
head(all_arg_sum_rpkm_HKH_CNH_CNS)

install.packages('ggtern')
library(ggtern)
fig4b_p<-ggtern(data=all_arg_sum_rpkm_HKH_CNH_CNS, aes(x=ChinaHuman,y=ChinaSwine, z=HKHuman,color=ACT)) + 
  geom_point(aes(size=scale(HKHuman)))+
  #tern_limits(0.4,0.4,0.5)+
  theme_bw()+
  labs(size = 'Scaled RPKM in HK Human') +
  theme(legend.text=element_text(size=18,family="Times",face = "bold"),
        legend.title = element_text(size=18,family="Times",face = "bold"),
        axis.text.x = element_text(face='bold',family="Times",size = 16),
        axis.text.y=element_text(face='bold',family="Times",size = 16),
        axis.title = element_text(family="Times",face="bold",size = 18))+
  #ggtitle('Feldspar - Elkins and Grove 1990')+
  #annotate(geom = 'text',xmax = 40,ymax = 40,zmax = 50)+
  scale_color_manual(values = pal_simpsons("springfield")(16))+
  xlab("CNH") +                       #replace default axis labels
  ylab("CNS") +
  zlab("HKH")#+

fig4b_p

cn_hk_shared_args<-data.frame(intersect(intersect(colnames(HKHuman_all_arg_res_rpkm_clean)[5:314],
          colnames(ChinaHuman_all_arg_res_rpkm_clean)[5:224]),
          colnames(ChinaSwine_all_arg_res_rpkm_clean)[5:329]))

colnames(cn_hk_shared_args)<-"arg"
head(cn_hk_shared_args)
arg_ACT<-read.csv('arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
arg_ACT<-rownames_to_column(arg_ACT, var = "arg")
colnames(arg_ACT)<-c("arg",'ACT')
cn_hk_shared_args<-merge(arg_ACT,cn_hk_shared_args,'arg')
head(cn_hk_shared_args)
dim(cn_hk_shared_args)
HKH_CNS_shared_arg<-intersect(colnames(HKHuman_all_arg_res_rpkm_clean)[5:314],
                              colnames(ChinaSwine_all_arg_res_rpkm_clean)[5:329])
cn_hk_shared_args<-intersect(intersect(colnames(HKHuman_all_arg_res_rpkm_clean)[5:314],
                                                  colnames(ChinaHuman_all_arg_res_rpkm_clean)[5:224]),
                                        colnames(ChinaSwine_all_arg_res_rpkm_clean)[5:329])
#p1<-ggplot(cn_hk_shared_args, aes(ACT)) + 
#  geom_bar(stat="identity",color="black", fill=pal_npg(c("nrc"))(10)[2])+
#  theme_bw()+
#  scale_fill_manual(values = pal_npg(c("nrc"))(10))
#p1

df<-as.data.frame(table(cn_hk_shared_args$ACT))
colnames(df)<-c("ATC","count")
df<-df[df$count!="0",]
df<-df[order(-df$count),]
fig4c_p1<-ggplot(df, aes(x = reorder(ATC, -count), y = count))+ 
  geom_bar(stat="identity", fill=pal_npg(c("nrc"))(10)[4],alpha=0.6)+
  xlab("") + ylab("Number of ARGs")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,size = 16,face = 'bold',family = "Times"),
        axis.text.y=element_text(size = 16,face = 'bold',family = "Times"),
        axis.title.y = element_text(size = 18,face = 'bold',family = "Times"))
fig4c_p1
#China_HK_IS<-read.csv('China_HK_IS.txt',header = F,sep='\t',row.names = 1)
CN_HK_IS <- read.csv('./China_HK_IS.txt',sep="\t", header = T,row.names = 1,stringsAsFactors = FALSE)
CN_HK_IS_0.005<-CN_HK_IS %>%
  rownames_to_column('IS') %>%
  filter_if(is.numeric,all_vars(.>=0.005)) %>%
  column_to_rownames('IS')

fig4d_p3<-pheatmap(as.matrix(CN_HK_IS_0.005),fontsize = 18,fontfamily = "Times",
         color=(colorRampPalette(brewer.pal(9, "Reds"))(1000)),
         silent = FALSE)

ratio_tax_CNS_CNH_HKH <- read.table('ratio_tax_CNS_CNH_HKH.txt',sep="\t", header = T,
                                  stringsAsFactors = FALSE)
CN_HK_tax_IS_m <- melt(ratio_tax_CNS_CNH_HKH, id.vars = "Tax", variable.name="Host", value.name = "Size")

head(CN_HK_tax_IS_m)

fig4e_p4<-ggplot(CN_HK_tax_IS_m, aes(x = Host, y = Tax,col=Host)) +
  geom_point(aes(size = Size*40),alpha=0.6) + 
  scale_size(name="Ratio of taxon",
             range = range(CN_HK_tax_IS_m$Size*40),
             #breaks = c(0, , 19, 21),
             labels = c(0, 0.1, 0.2, 0.3,0.4)) +
  xlab("")+ylab("Genus")+
  theme_bw()+theme(axis.text.x = element_text(colour = "grey20", size = 18,face = "bold",family="Times",
                                              angle = 90, hjust = 0.5, vjust = 0.5),
                   axis.text.y = element_text(colour = "grey20",face = "bold.italic", size = 16,family="Times"),
                   text = element_text(size = 18,family="Times"))+
  scale_color_manual(values = pal_npg(c("nrc"))(10)[c(1,2,4)])
fig4e_p4

### network for CNS,CNH and HKH
suppressMessages(library(igraph))

head(arg_tax_is_plasmids_1174_fil)
dim(arg_tax_is_plasmids_1174_fil)
arg_tax_is_plasmids_1174_IS<-subset(arg_tax_is_plasmids_1174_fil,V10!='')
head(arg_tax_is_plasmids_1174_IS)
levels(all_arg_res_rpkm_with_region$region)
regions<-c("ChinaSwine","ChinaHuman","HKHuman")
arg_is_CN_HK<-subset(arg_tax_is_plasmids_1174_IS,region %in% regions) %>%
  select(c('arg', 'region',rpkm='V6',IS='V10')) %>%
  unite(arg_type, c(arg, region), sep=':') #,rpkm
head(arg_is_CN_HK)
dim(arg_is_CN_HK)

#s <- strsplit(arg_is_1174_clean$IS, split = ",") #%>%unlist() %>%unnest()
s<-lapply(arg_is_CN_HK$IS, function (x) unique(unlist(strsplit(as.character(x), ","))))
s
length(s)
arg_is_CN_HK_m<-data.frame(arg_type = rep(arg_is_CN_HK$arg_type, sapply(s, length)), IS = unlist(s))
dim(arg_is_CN_HK_m)
head(arg_is_CN_HK_m)

ar.count <- arg_is_CN_HK_m%>%
  count(arg_type,IS)%>%
  data.frame()
dim(ar.count)
head(ar.count)

cluster.count<-arg_is_CN_HK_m%>%
  group_by(arg_type) %>%   ## group by annotation
  mutate(cluster.size=length(arg_type)) %>% 
  ungroup()%>%
  data.frame()
head(cluster.count)

regions<-c("ChinaSwine","ChinaHuman","HKHuman")
arg_is_CN_HK_rpkm<-subset(arg_tax_is_plasmids_1174_IS,region %in% regions) %>%
  select(c('arg', 'region',rpkm='V6',IS='V10')) %>%
  unite(arg_type, c(arg, region,rpkm), sep=':') #,rpkm
head(arg_is_CN_HK_rpkm)
dim(arg_is_CN_HK_rpkm)

#s <- strsplit(arg_is_1174_clean$IS, split = ",") #%>%unlist() %>%unnest()
s<-lapply(arg_is_CN_HK_rpkm$IS, function (x) unique(unlist(strsplit(as.character(x), ","))))
s
length(s)
arg_is_CN_HK_rpkm_m<-data.frame(arg_type = rep(arg_is_CN_HK_rpkm$arg_type, sapply(s, length)), IS = unlist(s))
dim(arg_is_CN_HK_rpkm_m)
head(arg_is_CN_HK_rpkm_m)

arg_rpkm<-arg_is_CN_HK_rpkm_m%>%
  mutate(rpkm=str_split_fixed(arg_type, ":", 3)[,3])%>%
  mutate(arg_type=str_remove(arg_type, ':[0-9].*+$'))%>%
  group_by(arg_type, IS) %>%
  summarise(arg_size=sum(as.numeric(rpkm)))%>%
  data.frame()
head(data.frame(arg_rpkm))
dim(arg_rpkm)

IS_rpkm<-arg_is_CN_HK_rpkm_m%>%
  mutate(rpkm=str_split_fixed(arg_type, ":", 3)[,3])%>%
  mutate(arg_type=str_remove(arg_type, ':[0-9].*+$'))%>%
  group_by(IS) %>%
  summarise(IS_size=sum(as.numeric(rpkm)))%>%
  data.frame()
head(IS_rpkm)
dim(IS_rpkm)

arg_is_CN_HK_dat<-merge(ar.count,cluster.count,by=c(1,2))%>%
  merge(arg_rpkm,by=c(1,2))%>%
  filter(cluster.size>=5 ) %>%
  mutate(score=n/cluster.size) %>%
  filter(score>= 0.6 )
dim(arg_is_CN_HK_dat)
head(arg_is_CN_HK_dat)
arg_is_CN_HK_dist<-distinct(arg_is_CN_HK_dat)
dim(arg_is_CN_HK_dist)
head(arg_is_CN_HK_dist)
write.table(arg_is_CN_HK_dist,'./arg_is_CN_HK_dist.txt', quote=F, sep='\t', row.names = F, col.names = T)
write.table(arg_is_CN_HK_dist,'./arg_is_CN_HK_cytoscape_dist.tsv', quote=F, sep='\t', row.names = F, col.names = T)

arg_ACT<-read.csv('arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
library(tibble)
arg_ACT<-rownames_to_column(arg_ACT, var = "arg")
colnames(arg_ACT)<-c("arg",'ACT')

arg_is_CN_HK_dist_ACT<-arg_is_CN_HK_dist%>%
  mutate(arg=str_split_fixed(arg_type, ":", 2)[,1])%>%
  merge(arg_ACT,'arg')%>%
  select(arg_type,IS,n,cluster.size,score,ACT,arg_size)%>%
  unite(arg_type,arg_type,ACT,arg_size,sep=':')
head(arg_is_CN_HK_dist_ACT)

g <- graph_from_data_frame(arg_is_CN_HK_dist_ACT, directed=FALSE)
V(g)$type <- 'NA'
V(g)$type[!(str_detect(V(g)$name, ':'))] <- 'IS'
V(g)$type[str_detect(V(g)$name, 'ChinaHuman')] <- 'ChinaHuman'
V(g)$type[str_detect(V(g)$name, 'ChinaSwine')] <- 'ChinaSwine'
V(g)$type[str_detect(V(g)$name, 'HKHuman')] <- 'HKHuman'

V(g)$lab <- V(g)$name
V(g)$color <- ifelse(V(g)$type=='IS', 'IS', ((str_split_fixed((V(g)$name),":", 4)[,3])))
V(g)$color
V(g)$lab <- ifelse((str_detect(V(g)$name, ':')),((str_split_fixed((V(g)$name),":", 2)[,1])),(V(g)$name))
#V(g)$type <- str_detect(V(g)$name, 'Human')
#V(g)$color <- ifelse(V(g)$type, 'ARGs', 'IS')
V(g)$lab[str_detect(V(g)$lab, 'ISCco2')] <- 'Tn5541/Tn5543'
V(g)$lab
E(g)$weight <- arg_is_CN_HK_dist_ACT$score
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
  scale_shape_manual(values=c(15,16,17,8),name='Host')+
  theme(legend.text=element_text(size=16,family="Times",face = "bold"),
        legend.title = element_text(size=16,family="Times","bold"))


cowplot::plot_grid(venn.plot,p1,p3$gtable,p4,labels = c('a.', 'b.','c.','d.'),label_size = 20,rel_heights = c(4,6))

top_row <- plot_grid(venn.plot,p1, labels = c('a.', 'b.'), label_size = 20)
bottom_row<-plot_grid(p3$gtable,p4, labels = c('c.', 'd.'), label_size = 20,rel_widths=c(12,8))
plot_grid(top_row, bottom_row, label_size = 12, ncol = 1,
                    rel_heights = c(1.5,2))
#left_col
plot_grid(left_col, p3$gtable,labels = c('','d.'), label_size = 12, nrow = 1,rel_widths=c(16,4))

