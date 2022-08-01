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
library(tibble)
setwd("~/Desktop/resistome")

#Gene Rarefaction
#generate gene rarefaction plots
#Load generic libraries
# https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
source('configuration.r')
#Load specific libraries

library(matrixStats)
library(iNEXT)
library(betapart)
library(tibble)
library(ggsci)
library(cowplot)

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

###subset samples referring to host
all_arg_res_rpkm_human<-subset(all_arg_res_rpkm_with_region,Host == "Human")
all_arg_res_rpkm_human<-all_arg_res_rpkm_human[,-c(2,3,4,868,869)]
rownames(all_arg_res_rpkm_human)<-all_arg_res_rpkm_human[,1]
all_arg_res_rpkm_human<-all_arg_res_rpkm_human[,-1]
class(all_arg_res_rpkm_human)
dim(all_arg_res_rpkm_human)
head(all_arg_res_rpkm_human)
all_arg_res_rpkm_human_clean <- ifelse(all_arg_res_rpkm_human>0,1,0) 
dim(all_arg_res_rpkm_human_clean)
head(all_arg_res_rpkm_human_clean)
all_arg_res_rpkm_human_clean_t<-t(all_arg_res_rpkm_human_clean)
class(all_arg_res_rpkm_human_clean_t)
dim(all_arg_res_rpkm_human_clean_t)
head(all_arg_res_rpkm_human_clean_t)
#all_arg_res_rpkm_human_clean_t<-data.frame(all_arg_res_rpkm_human_clean_t,check.names = FALSE)
typeof(all_arg_res_rpkm_human_clean_t)
all_arg_res_rpkm_human_clean_t

all_arg_res_rpkm_swine<-subset(all_arg_res_rpkm_with_region,Host == "Swine")
all_arg_res_rpkm_swine<-all_arg_res_rpkm_swine[,-c(2,3,4,868,869)]
rownames(all_arg_res_rpkm_swine)<-all_arg_res_rpkm_swine[,1]
all_arg_res_rpkm_swine<-all_arg_res_rpkm_swine[,-c(1)]
all_arg_res_rpkm_swine_clean <- ifelse(all_arg_res_rpkm_swine>0,1,0) 
all_arg_res_rpkm_swine_clean_t<-t(all_arg_res_rpkm_swine_clean)
#all_arg_res_rpkm_swine_mat[all_arg_res_rpkm_swine_mat !=0]<-1
dim(all_arg_res_rpkm_swine)
head(all_arg_res_rpkm_swine)

all_arg_res_rpkm_chicken<-subset(all_arg_res_rpkm_with_region,Host == "Chicken")
all_arg_res_rpkm_chicken<-all_arg_res_rpkm_chicken[,-c(2,3,4,868,869)]
rownames(all_arg_res_rpkm_chicken)<-all_arg_res_rpkm_chicken[,1]
all_arg_res_rpkm_chicken<-all_arg_res_rpkm_chicken[,-1]
head(all_arg_res_rpkm_chicken)
all_arg_res_rpkm_chicken_clean <- ifelse(all_arg_res_rpkm_chicken>0,1,0) 
all_arg_res_rpkm_chicken_clean_t<-t(all_arg_res_rpkm_chicken_clean)
dim(all_arg_res_rpkm_chicken_clean_t)

all_arg_res_rpkm_clean_m_list<-list(Human=all_arg_res_rpkm_human_clean_t,
                                          Swine=all_arg_res_rpkm_swine_clean_t,
                                          Chicken=all_arg_res_rpkm_chicken_clean_t)
all_out_all <- iNEXT(all_arg_res_rpkm_clean_m_list,datatype="incidence_raw", endpoint = 2000)
all_out_all
p2<-ggiNEXT(all_out_all)+
  labs(x="Number of samples", y="Species diversity") +
  scale_color_manual(values=pal_npg(c("nrc"))(10)[c(1,4,3)])+
  theme_bw()+
  theme(legend.position = "", ###legend.position = "bottom",
        legend.title=element_blank(),
        text=element_text(size=18,family = "Times"),
        legend.box = "vertical") 

all_arg_res_rpkm_human_clean_m_list

df_human<-ggplot2::fortify(all_out_human,type=1)
df_human[,"site"]<-"human"
ggiNEXT(all_out_human)

df_swine<-fortify(all_out_swine,type=1)
df_swine[,"site"]<-"swine"
ggiNEXT(all_out_swine)

all_arg_res_rpkm_chicken_mat
all_out_chicken <- iNEXT(all_arg_res_rpkm_chicken_mat, datatype="incidence_raw", endpoint = 2000)
df_swine<-fortify(all_out_chicken,type=1)
df_chicken[,"site"]<-"chicken"
ggiNEXT(all_out_chicken)

df_swine_human<-rbind(df_swine,df_human)
df_swine_human.point <- df_swine_human[which(df_swine_human$method=="observed"),]
df_swine_human.line <- df_swine_human[which(df_swine_human$method!="observed"),]
df_swine_human.line$method <- factor(df_swine_human.line$method, 
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))
df_swine					 

p2<-ggplot(df_swine_human, aes(x=x, y=y, colour=site)) + 
  geom_point(aes(shape=site), size=5, data=df_swine_human.point) +
  geom_line(aes(linetype=method), lwd=1.5, data=df_swine_human.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=site, colour=NULL), alpha=0.2) +
  labs(x="Number of samples", y="Species diversity") +
  scale_color_manual(values=pal_npg(c("nrc"))(10)[c(1,4)])+
  theme_bw()+
  theme(legend.position = "", ###legend.position = "bottom",
        legend.title=element_blank(),
        text=element_text(size=18,family = "Times"),
        legend.box = "vertical") 
p2
df_swine_human

install.packages("betapart")
library(betapart)
library(tibble)
library(ggsci)

par(mfrow=c(2,2))
##read data
all_arg_res_rpkm <- read.csv("../../microbiome/data/resfinder_all/1174all_arg_res_rpkm.txt",sep="\t", 
                             header = T, row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
all_arg_res_rpkm_t<-t(all_arg_res_rpkm)
all_arg_res_rpkm_t<-data.frame(all_arg_res_rpkm_t,check.names = FALSE)
all_arg_res_rpkm_t<-rownames_to_column(all_arg_res_rpkm_t, var = "sample") #### tibble is needed
sample_region<-read.csv('1174_sample.txt',header = F,sep='\t',row.names = 1,check.names = FALSE)
sample_region<-rownames_to_column(sample_region, var = "sample")
colnames(sample_region)<-c("sample",'region','type','nation')
all_arg_res_rpkm_with_region<-merge(sample_region,all_arg_res_rpkm_t,'sample')

all_arg_res_rpkm_human<-subset(all_arg_res_rpkm_with_region,type == "Human")
all_arg_res_rpkm_human<-all_arg_res_rpkm_human[,-c(2,3,4)]
rownames(all_arg_res_rpkm_human)<-all_arg_res_rpkm_human[,1]
all_arg_res_rpkm_human<-all_arg_res_rpkm_human[,-1]
all_arg_res_rpkm_human<-t(all_arg_res_rpkm_human)
all_arg_res_rpkm_human_mat<-all_arg_res_rpkm_human
all_arg_res_rpkm_human_mat[all_arg_res_rpkm_human_mat !=0]<-1

all_arg_res_rpkm_swine<-subset(all_arg_res_rpkm_with_region,type == "Swine")
all_arg_res_rpkm_swine<-all_arg_res_rpkm_swine[,-c(2,3,4)]
rownames(all_arg_res_rpkm_swine)<-all_arg_res_rpkm_swine[,1]
all_arg_res_rpkm_swine<-all_arg_res_rpkm_swine[,-1]
all_arg_res_rpkm_swine<-t(all_arg_res_rpkm_swine)
all_arg_res_rpkm_swine_mat<-all_arg_res_rpkm_swine
#all_arg_res_rpkm_swine_mat<-data.frame(all_arg_res_rpkm_swine_mat,check.names = FALSE)
all_arg_res_rpkm_swine_mat[all_arg_res_rpkm_swine_mat !=0]<-1

# Resample 100 times the multiple-site dissimilarities
# for 10 countries.
beta_human<-beta.sample(t(all_arg_res_rpkm_human_mat), index.family="sorensen",sites=10, samples=100)
beta_swine<-beta.sample(all_arg_res_rpkm_swine_mat, index.family="sorensen",sites=10, samples=100)
beta_chicken<-beta.sample(all_arg_res_rpkm_chicken_mat, index.family="sorensen",sites=10, samples=100)
# Plot the distributions of beta.SIM in Southern Europe (red) 
# and Northern Europe (blue)
cols=pal_simpsons("springfield")(16)
cols=pal_npg(c("nrc"))(10)
#show_col(cols)
plot(density(beta_human$sampled.values$beta.SIM), col=cols[1], xlim=c(0,1),lwd=3,main="Distributions of beta.SIM")
lines(density(beta_swine$sampled.values$beta.SIM), col=cols[2],lwd=3)
legend(1, 95, legend=c("Human", "Swine"),col=cols[15,16], cex=0.8)

human_beta<-data.frame(beta_human$sampled.values)%>%
  mutate(host="human")
swine_beta<-data.frame(beta_swine$sampled.values)%>%
  mutate(host="swine")
chicken_beta<-data.frame(beta_chicken$sampled.values)%>%
  mutate(host="chicken")

all_swine_human_beta<-rbind(human_beta,swine_beta,chicken_beta)

p3<-ggplot(all_swine_human_beta,aes(x=beta.SIM,fill=host)) +
  geom_density(alpha=0.6)+
  scale_fill_manual(values = pal_npg('nrc')(10)[c(1,4,3)])+
  theme_bw()+
  labs(fill="Host")+
  theme(axis.text.y = element_text(size = 12,family = "Times"),axis.text.x = element_text(size = 12,family = "Times"),
                         axis.title.x = element_text(size = 16,family = "Times"),axis.title.y = element_text(size = 16,family = "Times"),
                         legend.text=element_text(size=14,family = "Times"),legend.title=element_text(size=16,family = "Times"),
        legend.position = "")
p3
p4<-ggplot(all_swine_human_beta,aes(x=beta.SNE,fill=host)) +
  geom_density(alpha=0.6)+
  scale_fill_manual(values = pal_npg('nrc')(10)[c(1,4,3)])+
  theme_bw()+
  labs(fill="Host")+
  theme(axis.text.y = element_text(size = 12,family = "Times"),axis.text.x = element_text(size = 12,family = "Times"),
        axis.title.x = element_text(size = 16,family = "Times"),axis.title.y = element_text(size = 16,family = "Times"),
        legend.text=element_text(size=14,family = "Times"),legend.title=element_text(size=16,family = "Times"),
        legend.position = "")
p4
p5<-ggplot(all_swine_human_beta,aes(x=beta.SOR,fill=host)) +
  geom_density(alpha=0.6)+
  scale_fill_manual(values = pal_npg('nrc')(10)[c(1,4,3)])+
  theme_bw()+
  labs(fill="Host")+
  theme(axis.text.y = element_text(size = 12,family = "Times"),axis.text.x = element_text(size = 12,family = "Times"),
        axis.title.x = element_text(size = 16,family = "Times"),axis.title.y = element_text(size = 16,family = "Times"),
        legend.text=element_text(size=14,family = "Times"),legend.title=element_text(size=16,family = "Times"),
        legend.position = "")
p5

head(all_arg_res_rpkm_with_region)
dim(all_arg_res_rpkm_with_region)
#HK_all_arg_res_rpkm<-HK_all_arg_res_rpkm[, colSums(HK_all_arg_res_rpkm != 0) > 0]
all_arg_res_rpkm_with_region$arg_n<-rowSums(all_arg_res_rpkm_with_region[5:867] !=0)
#samples1174_ARG_count<-read.table('samples1174_ARG_count.txt', sep='\t',head=TRUE)
all_arg_res_rpkm_with_region$Host<-as.factor(all_arg_res_rpkm_with_region$Host)
p1<-ggplot(all_arg_res_rpkm_with_region, aes(arg_n,fill=Host)) + 
  geom_histogram(color="black",alpha=0.6,binwidth = 10)+
  labs(y="Number of samples", x='Number of ARGs')+
  scale_fill_manual(values = pal_npg('nrc')(10)[c(1,4,3)])+
  theme_bw()+
  labs(fill="Host")+
  theme(axis.text.y = element_text(size = 12,family = "Times"),
        axis.text.x = element_text(size = 12,family = "Times"),
        axis.title.x = element_text(size = 16,family = "Times"),
        axis.title.y = element_text(size = 16,family = "Times"),
        legend.text=element_text(size=14,family = "Times"),
        legend.title=element_text(size=16,family = "Times"))
p1

top_row <- plot_grid(p1,p2,labels = c('a.', 'b.'), label_size = 20,rel_widths=c(8,12))
bottom_row<-plot_grid(p3,p4,p5,labels = c('c.', 'd.','e.'), label_size = 20,nrow = 1)
plot_grid(top_row, bottom_row, label_size = 20, ncol = 1,
                    rel_heights = c(2,2))

ggsave("./microbiota_new_figures/beta-diversity.pdf", height = 10, width = 10)

# Compute the p-value of difference in beta.SIM between South and North 
# (i.e. the probability of finding in the North a higher value than 
# in the South)

p.value.beta.SIM<-length(which(beta_human$sampled.values$beta.SIM<
                                 beta_swine$sampled.values$beta.SIM))/100

p.value.beta.SIM
# The result is 0 and we used 100 samples, so p<0.01
