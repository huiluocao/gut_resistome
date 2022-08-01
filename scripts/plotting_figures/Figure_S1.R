#---------------------------------------------------------------------------------------------------------
library(maps)
library(ggplot2)
library(RColorBrewer)
library(ggforce)
library(reshape2)
library(dplyr)
library(grid)
library(ggsci)

world_map <- map_data("world")
unique(world_map['region'])
df<-read.csv("./data/sample_stat.csv",header=TRUE)
#df<-read.csv("../microbiome/old/microbiota_R/sample_stat.csv",header=TRUE)
Region_Center<-aggregate(cbind(long,lat)~region,world_map,mean)

df_join<-left_join(df,Region_Center,by='region')
df_join[df_join$region=='Hong Kong',c('lat','long')]=c(22.28552, 114.15769)


df<-melt(df_join,id.vars=c('region','long','lat'))
df
df$start<- rep(c(-pi/2, pi/2), each=nrow(df)/2)
df<-transform(df,plus=ifelse(variable=='Human',1,-1))
r <- 8
min_r<-0.5
scale <- r/max(sqrt(df$value))

library("scales")
show_col(pal_npg('nrc')(10))

p1_map<-ggplot() +
  geom_path(data=world_map,aes(x=long,y=lat,group=group),colour="gray80",size=.2)+ ###colour="gray80"
  geom_arc_bar(data=df,aes(x0 = long, y0 =lat, r0 = 0, r = sqrt(value)*scale+min_r,
                           start = start, end = start + pi, fill = variable),
              size=0.2,alpha=0.8) + ## color = "black",
  geom_text(data=df,aes(label = value,size=value, x = long, y = lat + plus*(scale*sqrt(value)/2*0.9+0.25)),
            color='black')+  
  scale_fill_manual(values = pal_npg('nrc')(10)[c(1,4,3)])+
  scale_x_continuous(limits=c(-10,30),expand=c(0,0))+
  scale_y_continuous(limits=c(30,70),expand=c(0,0))+#breaks=(-3:3)*30
  guides(fill = FALSE,size=FALSE) +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_void()+
  theme(panel.background =element_rect(fill=NA,color='black') )
p1_map


p2_map<-ggplot() + 
  geom_path(data=world_map,aes(x=long,y=lat,group=group),colour="gray80",size=.2)+
  geom_arc_bar(data=df,aes(x0 = long, y0 =lat, r0 = 0, r = sqrt(value)*scale+min_r,
                           start = start, end = start + pi, fill = variable),
               color = "black",size=0.2,alpha=0.8) +
  
  geom_text(data=df,aes(label = value,size=value, x = long, y = lat + plus*(scale*sqrt(value)/2*0.9+0.25)),
             color='black')+  
  scale_size(range=c(0,4))+#scale_x_continuous(limits=c(-10,30),expand=c(0,0))+
  scale_fill_manual(values = pal_npg('nrc')(10)[c(1,4)])+
  #scale_y_continuous(limits=c(30,70),expand=c(0,0))+#breaks=(-3:3)*30
  guides(fill = guide_legend(title = "Host"),size=FALSE) +
  xlab("Longitude") + 
  ylab("Latitude") +
  coord_fixed() +
  theme_test()+
  theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16),
        legend.text=element_text(size=14),legend.title=element_text(size=16))
p2_map

subvp<-viewport(x=0.4,0.5,width=0.3,height=0.4)
p2_map
print(p1_map,vp=subvp)
  
ggsave("./microbiota_new_figures/map.pdf", height = 10, width = 10)
