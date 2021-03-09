library(vegan)
library(ggpubr)
####shannon and simpson indices for each subject.
dim(all_arg_res_rpkm_with_region)
head(all_arg_res_rpkm_with_region)
all_arg_res_rpkm_with_region$shanno<-diversity(all_arg_res_rpkm_with_region[,5:867],index = "shannon", MARGIN=1,base=exp(1))
all_arg_res_rpkm_with_region$simpson<-diversity(all_arg_res_rpkm_with_region[,5:867],index = "simpson", MARGIN=1,base=exp(1))
s4_p<-ggboxplot(all_arg_res_rpkm_with_region, x="Type", y="shanno",x.text.angle=90,
              add = "jitter", color = "Type",remove = c('IcelandHuman'),font.label=c(14, "bold"),
              xlab = FALSE,ylab = "Shannon index",yticks.by=0.8,xtickslab = FALSE,ggtheme = theme_bw()) ### facet.by = "type",
s4_p<-ggpar(s4_p,xtickslab = FALSE,legend = "none",font.y = c(14, "bold"))
s4_p<-s4_p+rremove("x.text")
s5_p<-ggboxplot(all_arg_res_rpkm_with_region, x="Type", y="simpson",
              add = "jitter", color = "Type",remove = c('IcelandHuman'),ylim=c(0,1.2),
              ylab = "Simpson index",ggtheme = theme_bw()) ###facet.by = "type",
s5_p<-ggpar(s5_p,x.text.angle=90,legend='none',xlab = "",font.xtickslab=c(14, "bold"),font.y = c(14, "bold"))
ggarrange(s4_p,s5_p,nrow = 2,heights = 3:4)
