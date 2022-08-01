library(dplyr)
library(tibble)
library(ggplot2)
library(ggsci)
library(hrbrthemes)

sample_region<-read.csv('1174_sample.txt',header = F,sep='\t',row.names = 1,check.names = FALSE)
sample_region<-rownames_to_column(sample_region, var = "sample")
colnames(sample_region)<-c("sample",'region','type',"nation")
sample_region<-sample_region[,c(1,2)]
head(sample_region)
dim(sample_region)

arg_ACT<-read.csv('data/arg_ATC_new.txt',header = F,sep='\t',row.names = 1)
head(arg_ACT)
library(tibble)
arg_ACT<-rownames_to_column(arg_ACT, var = "arg")
colnames(arg_ACT)<-c("arg",'ACT')
arg_sample_t_class<-merge(arg_ACT,arg_sample_t,by.x = 'arg',by.y = 'row.names')
head(arg_sample_t_class)
dim(arg_sample_t_class)

library(plyr)
library(data.table)
samples1174_ARG_count<-data.frame(t(numcolwise(nonzero)(all_arg_res_rpkm_1487)))
colnames(samples1174_ARG_count)<-c("n")
head(samples1174_ARG_count)
dim(samples1174_ARG_count)

arg_tax_is_plasmids_1174<-read.table('../../microbiome/data/arg_tax_is_plasmids_1174.txt',
                                     stringsAsFactors = FALSE, sep='\t', head=F, quote = "#")
##arg_tax_is_plasmids_1174_IS<-subset(arg_tax_is_plasmids_1174,!is.na(arg_tax_is_plasmids_1174$V10))
dim(arg_tax_is_plasmids_1174)
head(arg_tax_is_plasmids_1174)
arg_tax_is_plasmids_1174_fil<-merge(sample_region,arg_tax_is_plasmids_1174,by.x='sample',by.y='V1')
arg_tax_is_plasmids_1174_fil<-merge(arg_ACT,arg_tax_is_plasmids_1174_fil,by.x='arg',by.y='V7')
head(arg_tax_is_plasmids_1174_fil)
dim(arg_tax_is_plasmids_1174_fil)

library(dplyr)
all_arg_sum_1174<-arg_tax_is_plasmids_1174_fil %>% 
  group_by(sample) %>% 
  summarise(all_arg = sum(V6))%>%
  data.frame()%>%
  mutate(log2all_arg=log2(all_arg))
dim(all_arg_sum_1174)
head(all_arg_sum_1174)
all_arg_sum_1174_with_is<-arg_tax_is_plasmids_1174_fil %>%
  subset(V10!='') %>%
  group_by(sample) %>% 
  summarise(all_arg_is = sum(V6))%>%
  data.frame()%>%
  mutate(log2all_arg_is=log2(all_arg_is))
head(all_arg_sum_1174_with_is)
dim(all_arg_sum_1174_with_is)
all_arg_sum_1174_fil<-all_arg_sum_1174[which(all_arg_sum_1174$sample %in% all_arg_sum_1174_with_is$sample),]
head(all_arg_sum_1174_fil)
dim(all_arg_sum_1174_fil)

all_arg_sum_1174_model<-merge(all_arg_sum_1174_fil,all_arg_sum_1174_with_is,by='sample') %>%
  mutate(log2ratio=log2all_arg_is/log2all_arg)%>%
  merge(sample_region,by='sample')%>%
  merge(samples1174_ARG_count,by.x = 'sample',by.y = 0)%>%
  filter(log2ratio>=0)

all_arg_sum_1174_model_scale<-merge(all_arg_sum_1174_fil,all_arg_sum_1174_with_is,by='sample') %>%
  mutate(scale_arg=scale(all_arg))%>%
  mutate(ratio=(all_arg_is/all_arg))%>%
  merge(sample_region,by='sample')%>%
  merge(samples1174_ARG_count,by.x = 'sample',by.y = 0)

fig5_p1<-ggstatsplot::grouped_ggscatterstats(
  data = all_arg_sum_1174_model_scale,
  x=n,
  y=ratio,
  grouping.var = Host,
#  size=scale_arg,color=nation,shape=type,
  type = "robust", # type of test that needs to be run
  conf.level = 0.99, # confidence level
  xlab = "Number of ARGs", # label for x axis
  ylab = "Ratio of ARGs associated with IS", # label for y axis
  label.var = "sample", # variable for labeling data points
  label.expression = "ratio > 0.8 & log2all_arg_is > 10", # expression that decides which points to label
  line.color = "yellow", # changing regression line color line
  #title = "",# title text for the plot
  title.prefix = "Host",
#  caption = expression( # caption text for the plot
#    paste(italic("type"), ": ARGs")
#  ),
  ggtheme = ggthemes::theme_few(), # choosing a different theme,#hrbrthemes::theme_ipsum_ps()
  ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
  marginal.type = "density", # type of marginal distribution to be displayed
  xfill = "#0072B2", # color fill for x-axis marginal distribution
  yfill = "#009E73", # color fill for y-axis marginal distribution
  xalpha = 0.6, # transparency for x-axis marginal distribution
  yalpha = 0.6, # transparency for y-axis marginal distribution
  centrality.para = "median", # central tendency lines to be displayed
  messages = FALSE # turn off messages and notes
#plotgrid.args = list(nrow = 2)
) +ggplot2::theme(plot.title = element_text(size = 14, face = "bold",family = "Times"),
                 text = element_text(size = 14,family = "Times"),
                 axis.title = element_text(face="bold",family = "Times"),
                 axis.title.x = element_text(face="bold",family = "Times"),
                 axis.text.x = element_text(size = 14,family = "Times"))
#                 plotgrid.args = list(ncol = 3))
#                 legend.position = "none")
fig5_p1
