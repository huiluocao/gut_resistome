####some statisctics
library(plyr)
library(ggridges)
library(ggplot2)
library(gplots)
library(histogram)
library(tibble)
library(ggsci)

colors=pal_npg('nrc')(10)
pal_simpsons("springfield")(16)

dim(all_arg_res_rpkm_1487)
head(all_arg_res_rpkm_1487)
all_arg_res_rpkm_1487[is.na(all_arg_res_rpkm_1487)] <- 0
par(mfrow=c(1,3))
colSums(all_arg_res_rpkm_1487 !=0)
nonzero <- function(x) sum(x != 0)
numcolwise(nonzero)(all_arg_res_rpkm_1487)

#histogram(t(numcolwise(nonzero)(all_arg_res_rpkm_1487)),col=colors[1:11],
#          main="Frequency of ARGs",xlab = "Number of ARGs",ylab="Percent of samples") ### library(histogram),library(histogram)

samples1174_ARG_count<-data.frame(t(numcolwise(nonzero)(all_arg_res_rpkm_1487)))
colnames(samples1174_ARG_count)<-c("ARG_count")
head(samples1174_ARG_count)
write.table(samples1174_ARG_count,'./data/samples1487_ARG_count.txt', sep='\t', quote=F, row.names = T, col.names = T)
sample_region1<-read.csv('./data/1487_sample1.txt',header = T,sep='\t',check.names = FALSE,stringsAsFactors=T)
head(sample_region1)
samples1174_ARG_count<-merge(sample_region1,samples1174_ARG_count, by="row.names")
colnames(samples1174_ARG_count)<-c("sample",'Type','Host',"Region",'Gb','ARGs_per_Gbs','ARG_count')
#samples1174_ARG_count<-samples1174_ARG_count%>%
#  mutate(ARGs_per_Gbs=as.numeric(ARG_count)/as.numeric(Gb))
  
figS2_p <- ggplot(data =sample_region1, aes(y=Host,x=ARG_per_Gb, fill=Host, height = ..density..)) + 
  #geom_density_ridges(trim = TRUE, alpha=0.7) +
  stat_density_ridges(alpha = 0.3,  scale=0.98) +
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE) +
  scale_fill_npg(guide=FALSE)+
  theme_bw()+
  theme(axis.text.y = element_text(face='bold',size = 18))  + 
  scale_x_log10() + 
  labs(y=NULL, x='Number of ARGs per Gb reads') +
  coord_cartesian(ylim =c(1.5,3.5))
figS2_p
