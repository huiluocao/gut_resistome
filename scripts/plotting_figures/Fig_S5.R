library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(stringr) ### to use str_replace
library(tidyr) ## to use spread
library(tidyverse)
library(reshape2)

#Load functions
### for all ARGs in all 1487 samples
### for all ARGs in all 697 human samples
### for all ARGs in all 313 chicken samples
all_arg_res_rpkm_1174 <- read.csv("./data/1174all_arg_res_rpkm.txt",sep="\t", 
                                  header = T, row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
head(all_arg_res_rpkm_1174)
dim(all_arg_res_rpkm_1174)
row.names(all_arg_res_rpkm_1174)
all_arg_res_rpkm_313 <- read.csv("./data/all_chicken_arg_res_rpkm.txt",sep="\t", 
                                 header = T, row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
head(all_arg_res_rpkm_313)
dim(all_arg_res_rpkm_313)

all_arg_res_rpkm_1487 <- merge(all_arg_res_rpkm_1174,all_arg_res_rpkm_313,by=0,all=TRUE)
head(all_arg_res_rpkm_1487)
dim(all_arg_res_rpkm_1487)
all_arg_res_rpkm_1487<-column_to_rownames(all_arg_res_rpkm_1487,var = 'Row.names')
#full_join(df1, df2, by = "wpt")

all_arg_res_rpkm_1487_t<-t(all_arg_res_rpkm_1487)
head(all_arg_res_rpkm_1487_t)
dim(all_arg_res_rpkm_1487_t)
all_arg_res_rpkm_1487_t[is.na(all_arg_res_rpkm_1487_t)] <- 0
all_arg_res_rpkm_1487_t<-data.frame(all_arg_res_rpkm_1487_t,check.names = FALSE)
#disttransform(all_arg_res_rpkm, method="hellinger")
dim(all_arg_res_rpkm_1487_t)

all_arg_res_rpkm_1487_t<-rownames_to_column(all_arg_res_rpkm_1487_t, var = "sample") #### tibble is needed

sample_region<-read.csv('./data/1487_sample.txt',header = F,sep='\t',row.names = 1,check.names = FALSE,stringsAsFactors=T)
sample_region<-rownames_to_column(sample_region, var = "sample")
colnames(sample_region)<-c("sample",'Type','Host',"Region")
head(sample_region)
dim(sample_region)
#write.csv(sample_region,file = "sample_region.csv",row.names = TRUE, col.names = TRUE,quote=FALSE)

##merge ARGs abundance data with metadata
all_arg_res_rpkm_with_region<-merge(sample_region,all_arg_res_rpkm_1487_t,'sample')
dim(all_arg_res_rpkm_with_region)
head(all_arg_res_rpkm_with_region)
all_arg_res_rpkm_with_region$sum<-rowSums(all_arg_res_rpkm_with_region[,5:867])
all_arg_res_rpkm_with_region$log2sum<-log2(all_arg_res_rpkm_with_region$sum)
levels(all_arg_res_rpkm_with_region$Region)

##Merge data
dim(all_arg_res_rpkm_with_region)
head(all_arg_res_rpkm_with_region)
levels(all_arg_res_rpkm_with_region$Type)

arg_data<-all_arg_res_rpkm_with_region[,-c(1,3,4)]
arg_data_type<-aggregate(.~Type,arg_data, sum)
head(arg_data_type)
dim(arg_data_type)
arg_data_type$Type
rownames(arg_data_type)<-arg_data_type$Type
arg_data_type<-arg_data_type[,-1]
arg_data_type_t<-data.frame(t(arg_data_type),check.names = FALSE)
head(arg_data_type_t)
dim(arg_data_type_t)
arg_data_type_t<-arg_data_type_t[!(row.names(arg_data_type_t) %in% c('sum')), ]

arg_data_type_t<-rownames_to_column(arg_data_type_t, var = "Anti")
arg_data_type_melted<-melt(arg_data_type_t,id=c("Anti"), value.name="n")
head(arg_data_type_melted)
dim(arg_data_type_melted)

arg_ACT<-read.csv('./data/arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
arg_ACT<-rownames_to_column(arg_ACT, var = "Anti")
colnames(arg_ACT)<-c("Anti",'Anti_type')
colnames(arg_data_type_melted)<-c('Anti',"Sample_type","n")
arg_data_type_ATC<-merge(arg_ACT,arg_data_type_melted,'Anti')
arg_data_type_ATC_log2n<-mutate(arg_data_type_ATC,log2n=log2(arg_data_type_ATC$n))
arg_data_type_ATC_log2n$log2n[is.infinite(arg_data_type_ATC_log2n$log2n)]<-0
arg_data_type_ATC_log2n<-subset(arg_data_type_ATC_log2n,Sample_type!="IcelandHuman")
data.frame(arg_data_type_ATC_log2n)
arg_data_type_ATC_log2n$Anti_type<-as.factor(arg_data_type_ATC_log2n$Anti_type)
dim(arg_data_type_ATC_log2n)
head(arg_data_type_ATC_log2n)
unique(arg_data_type_ATC_log2n$Sample_type)

###Beta lactamase:
all_beta_lactamase<- subset(arg_data_type_ATC_log2n,Anti_type=="Beta-Lactams")
bla_class<-read.table('./data/Beta_lactamase_class.txt',header = F,sep='\t')
colnames(bla_class)<-c("Anti","bla_class")
all_beta_lactamase_class<-merge(bla_class,all_beta_lactamase,'Anti')

all_beta_lactamase_class.w <- subset(all_beta_lactamase_class,bla_class=="Class_C")%>%
  select(-Anti_type,-n,-bla_class)%>%
  spread(Sample_type, log2n, fill=0)%>%
  column_to_rownames("Anti")
dim(all_beta_lactamase_class.w)
head(all_beta_lactamase_class.w)
pheatmap(all_beta_lactamase_class.w,color=(colorRampPalette(brewer.pal(9, "Reds"))(1000)), main="Beta-lactamase Class C",
         legend_breaks = c(0,1, 2, 5, 10, 20),
         legend_labels = as.character(c(0, 1, 2, 5, 10, 20)),
         fontsize = 12, fontsize_row = 10, fontsize_col = 10,
         border_color="black",cluster_rows = TRUE, cluster_cols = TRUE, silent=FALSE)

all_beta_lactamase_class.w <- subset(all_beta_lactamase_class,bla_class=="Class_A")%>%
  select(-Anti_type,-n,-bla_class)%>%
  spread(Sample_type, log2n, fill=0)%>%
  column_to_rownames("Anti")
pheatmap(all_beta_lactamase_class.w,color=(colorRampPalette(brewer.pal(9, "Greens"))(1000)), main="Beta-lactamase Class A",
         legend_breaks = c(0,1, 2, 5, 10, 20),
         legend_labels = as.character(c(0, 1, 2, 5, 10, 20)),
         fontsize = 12, fontsize_row = 10, fontsize_col = 10,
         border_color="black",cluster_rows = TRUE, cluster_cols = TRUE, silent=FALSE)

all_beta_lactamase_class.w <- subset(all_beta_lactamase_class,bla_class=="Class_B1" | bla_class=="Class_B3")%>% ### no B2 was detected in all samples
  select(-Anti_type,-n,-bla_class)%>%
  spread(Sample_type, log2n, fill=0)%>%
  column_to_rownames("Anti")
pheatmap(all_beta_lactamase_class.w,color=(colorRampPalette(brewer.pal(9, "PuBu"))(1000)), main="Beta-lactamase Class B",
         legend_breaks = c(0,1, 2, 5, 10, 20),
         legend_labels = as.character(c(0, 1, 2, 5, 10, 20)),
         fontsize = 12, fontsize_row = 10, fontsize_col = 10,
         border_color="black",cluster_rows = TRUE, cluster_cols = TRUE, silent=FALSE)

all_beta_lactamase_class.w <- subset(all_beta_lactamase_class,bla_class=="Class_D")%>% 
  select(-Anti_type,-n,-bla_class)%>%
  spread(Sample_type, log2n, fill=0)%>%
  column_to_rownames("Anti")
pheatmap(all_beta_lactamase_class.w,color=(colorRampPalette(brewer.pal(9, "Purples"))(1000)), main="Beta-lactamase Class D",
         legend_breaks = c(0,1, 2, 5, 10, 20),
         legend_labels = as.character(c(0, 1, 2, 5, 10, 20)),
         fontsize = 12, fontsize_row = 10, fontsize_col = 10,
         border_color="black",cluster_rows = TRUE, cluster_cols = TRUE, silent=FALSE)


plot.heatmap <- function(anti_type, col='PuBu', cluster_rows=TRUE, cluster_cols=TRUE){
  arg_data_type_ATC_log2n.w <- subset(arg_data_type_ATC_log2n,Anti_type==anti_type) %>% 
    select(-Anti_type,-n) %>% 
    spread(Sample_type, log2n, fill=0) %>% 
    column_to_rownames("Anti")
  pheatmap(arg_data_type_ATC_log2n.w,color=(colorRampPalette(brewer.pal(9, col))(1000)), main=anti_type,
           legend_breaks = c(0,1, 2, 5, 10, 20),
           legend_labels = as.character(c(0, 1, 2, 5, 10, 20)),
           fontsize = 12, fontsize_row = 10, fontsize_col = 10,
           border_color="black",cluster_rows = cluster_rows, cluster_cols = cluster_cols, silent=FALSE)
}

#Fig 1e
plot.heatmap("Beta-lactamase", 'Reds') ### needs to be filtered
plot.heatmap("Aminoglycosides",'Blues')
plot.heatmap("Tetracyclines", 'Greens')
plot.heatmap("Colistins", 'Purples')
plot.heatmap("Fosfomycins","Reds")

unique(arg_data_type_ATC_log2n$Anti_type)
plot.heatmap("Fusidic acids","Reds") ### only ARGs in this class were detected in ItalySwine samples.
arg_data_type_ATC_log2n.w <- subset(arg_data_type_ATC_log2n,Anti_type=="Fusidic acid") %>% 
  select(-Anti_type,-n) %>% 
  spread(Sample_type, log2n, fill=0) %>% 
  column_to_rownames("Anti")
dim(arg_data_type_ATC_log2n.w)
head(arg_data_type_ATC_log2n.w)
pheatmap(arg_data_type_ATC_log2n.w,color=(colorRampPalette(brewer.pal(9, "Reds"))(1000)), main=anti_type,
         legend_breaks = c(0,1, 2, 5, 10, 20),
         legend_labels = as.character(c(0, 1, 2, 5, 10, 20)),
         fontsize = 12, fontsize_row = 10, fontsize_col = 10,
         border_color="black",cluster_rows = T, cluster_cols = T, silent=FALSE)

plot.heatmap("Rifampicins", 'RdPu')
plot.heatmap("Macrolides", 'YlGn')
plot.heatmap("Nitroimidazoles",'Blues')
plot.heatmap("Phenicols", 'Greys')
plot.heatmap("Quinolones", 'Oranges')
plot.heatmap("Sulphonamides", 'GnBu')
plot.heatmap("Trimethoprims", 'Greys')
plot.heatmap("Vancomycins", 'PuBuGn')

library(pdftools)
pdf_combine(c("./figures/suppls/Aminoglycoside_heatmap_all.pdf",
              "./figures/suppls/Beta_lactamase_class_A.pdf",
              "./figures/suppls/Beta_lactamase_class_B.pdf",
              "./figures/suppls/Beta_lactamase_class_C.pdf",
              "./figures/suppls/Beta_lactamase_class_D.pdf",
              "./figures/suppls/Colistin_heatmap_all.pdf",
              "./figures/suppls/Fosfomycin_heatmap_all.pdf",
              "./figures/suppls/Macrolide_heatmap_all.pdf",
              "./figures/suppls/Nitroimidazole_heatmap_all.pdf",
              "./figures/suppls/Phenicol_heatmap_all.pdf",
              "./figures/suppls/Quinolone_heatmap_all.pdf",
              "./figures/suppls/Rifampicin_heatmap_all.pdf",
              "./figures/suppls/Sulphonamide_heatmap_all.pdf",
              "./figures/suppls/Tetracycline_heatmap_all.pdf",
              "./figures/suppls/Trimethoprim_heatmap_all.pdf",
              "./figures/suppls/Vancomycins_heatmap_all.pdf"), 
            output = "./figures/Fig_S5.pdf")

### colors can be chosen here:
### https://www.rdocumentation.org/packages/RColorBrewer/versions/1.1-2/topics/RColorBrewer

#Suppl 2
s1 <- plot.heatmap("AGly")
s2 <- plot.heatmap("Tet", 'Greens')
s3 <- plot.heatmap("Phe", 'Purples')
s4 <- plot.heatmap("MLS", 'RdPu')
s5 <- plot.heatmap("Sul", 'YlGn')
s6 <- plot.heatmap("Gly",'Blues')
s7 <- plot.heatmap("Rif", 'Greys', cluster_rows = F, cluster_cols = F, log_mat = F)
s8 <- plot.heatmap("Tmt", 'Oranges')
s9 <- plot.heatmap("Flq", 'GnBu')
#Main

cowplot::plot_grid(m1$gtable)

plot.dat <- group_by(arg_data_region_ATCa, Sample_type) %>%
  tally() %>%
  merge(df.anti, by='Sample_type') %>%
  mutate(prev=n.y/n.x)

metadata <-read.table("./hospital_microbiome/metadata/illumina_metadata.txt",sep="\t",head=T) %>% 
  filter(Room_type != 'GIS')
anti <-read.table("../hospital_microbiome/tables/illumina_AR_gene_assignment.dat",sep="\t",head=T)
df.anti <- merge(metadata, anti, by.x='Site', by.y='Lib') %>% 
  mutate(Anti_type=gsub('.*_','',Anti)) %>%
  group_by(Sample_type,Anti_type,Anti,Library) %>%
  summarise() %>%
  group_by (Sample_type, Anti_type, Anti) %>%
  tally()

plot.dat <- group_by(metadata, Sample_type) %>%
  tally() %>%
  merge(df.anti, by='Sample_type') %>%
  mutate(prev=n.y/n.x) %>% 
  mutate(Sample_type= str_replace(Sample_type, "_"," ") %>% 
           str_replace("Door handle-interior", "Door Handle")) %>% 
  mutate(Anti= str_replace(Anti, "_[a-zA-Z]+$","") %>% 
           str_replace("_"," "))
#Helper function to plot

plot.heatmap <- function(anti_type, col='PuBu', cluster_rows=TRUE, cluster_cols=TRUE, log_mat=TRUE){
  plot.dat.w <- filter(plot.dat,Anti_type==anti_type) %>% 
    select(-Anti_type,-n.x,-n.y) %>% 
    spread(Sample_type, prev, fill=0) %>% 
    column_to_rownames("Anti") 
  if(!log_mat){
    pheatmap(plot.dat.w*100,color=(colorRampPalette(brewer.pal(9, col))(1000)), main=anti_type, 
             fontsize = 12, fontsize_row = 10, fontsize_col = 10,
             border_color="black",cluster_rows = cluster_rows, cluster_cols = cluster_cols, silent=TRUE)
  }else{
    pheatmap(log10(plot.dat.w*100+1),color=(colorRampPalette(brewer.pal(9, col))(1000)), main=anti_type,
             legend_breaks = log10(c(0,1, 2, 5, 10, 20, 50, 100, 120)+1),
             legend_labels = as.character(c(0, 1, 2, 5, 10, 20, 50, 100, 120)),
             fontsize = 12, fontsize_row = 10, fontsize_col = 10,
             border_color="black",cluster_rows = cluster_rows, cluster_cols = cluster_cols, silent=TRUE)
  }
}
###Run plot

#Fig 1e
m1 <- plot.heatmap("Bla", 'Reds', TRUE, TRUE)
#Suppl 2
s1 <- plot.heatmap("AGly")
s2 <- plot.heatmap("Tet", 'Greens')
s3 <- plot.heatmap("Phe", 'Purples')
s4 <- plot.heatmap("MLS", 'RdPu')
s5 <- plot.heatmap("Sul", 'YlGn')
s6 <- plot.heatmap("Gly",'Blues')
s7 <- plot.heatmap("Rif", 'Greys', cluster_rows = F, cluster_cols = F, log_mat = F)
s8 <- plot.heatmap("Tmt", 'Oranges')
s9 <- plot.heatmap("Flq", 'BuGn')
#Main

cowplot::plot_grid(m1$gtable)

ggsave('../plots/fig1e_antibiotics_profile.pdf', height = 7, width = 4)

#Supplementary figures

cowplot::plot_grid(
  cowplot::plot_grid(s1$gtable, s2$gtable, s3$gtable, s4$gtable, nrow=1),
  NULL,
  cowplot::plot_grid(s5$gtable, s6$gtable, s7$gtable, s8$gtable, s9$gtable, nrow=1, align='h'),
  NULL,
  ncol=1,
  rel_heights = c(4,0.5,2,0.5)
)
