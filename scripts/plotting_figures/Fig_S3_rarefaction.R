library(matrixStats)
library(iNEXT)
#Load functions
### for all ARGs in all 1487 samples
dim(all_arg_res_rpkm_1487)
head(all_arg_res_rpkm_1487)
all_arg_res_rpkm_1487[is.na(all_arg_res_rpkm_1487)] <- 0

all_arg_res_rpkm_1487_rarefac<-all_arg_res_rpkm_1487
all_arg_res_rpkm_1487_rarefac[all_arg_res_rpkm_1487_rarefac !=0]<-1
all_arg_res_rpkm_1487_rarefac
all_arg_res_rpkm_1487_rarefac_list<-list(All=as.matrix(all_arg_res_rpkm_1487_rarefac))
all_arg_res_rpkm_1487_rarefac_list            
all_out <- iNEXT(all_arg_res_rpkm_1487_rarefac_list, datatype="incidence_raw", endpoint = 5000)
g1 <- ggiNEXT(all_out, color.var="order") +
  geom_line(aes(y=all_out$AsyEst["Species Richness", "Estimator"]), size = 0.5, linetype = "dashed", color = "black") + 
  scale_y_continuous(breaks = c(seq(0 , 200, by=50), round(all_out$AsyEst["Species Richness", "Estimator"]))) +
  labs(title = "Rarefaction of all AR gene", x = "number of samples")
g1 <- g1 + guides(fill=FALSE)
g1 <- g1 + scale_color_discrete(labels = c("richness")) + 
  scale_shape_discrete(labels = c("richness"))
g1
index_out <- iNEXT(all_arg_res_rpkm_1487_rarefac_list, datatype="incidence_raw", q=c(1,2),endpoint = 5000)
g2 <- ggiNEXT(index_out, color.var="order") +
  geom_line(aes(y=index_out$AsyEst["Shannon diversity", "Estimator"]), size = 0.5, linetype = "dashed", color = "black") +
  geom_line(aes(y=index_out$AsyEst["Simpson diversity", "Estimator"]), size = 0.5, linetype = "dashed", color = "black") +
  scale_y_continuous(breaks = c(seq(0 , 100, by=10), round(index_out$AsyEst["Shannon diversity", "Estimator"]), round(index_out$AsyEst["Simpson diversity", "Estimator"])))+
  labs(title = "Rarefaction of all AR gene", y = "diversity index", x = "number of samples")
g2
g2 <- g2 + guides(fill=FALSE)
g2 <- g2 + scale_color_discrete(labels = c("shannon", "simpson")) +
  scale_shape_discrete(labels = c("shannon", "simpson"))
cowplot::plot_grid(g1, g2, nrow=2)
