#NON-METRIC MULTIDIMENSIONAL SCALING ANALYSIS
#FOR TETHYS DATASET AND THREE REGIONS (IRAN, SOUTH CHINA, PAKISTAN)
#RUN OCCURRENCE_ACQ.R FIRST
#MUST HAVE VEGAN PACKAGE DOWNLOADED

library(vegan)

#Runs NMDS with Bray-Curtis distance (default) but without row/column standardization (autotransform=F)
iran_MDS<-metaMDS(t(iran_ordered),autotransform=F) #performs non-metric multidimensional scaling
schina_MDS<-metaMDS(t(schina_ordered),autotransform=F)
pakistan_MDS<-metaMDS(t(pakistan_ordered),autotransform=F)
tethys_MDS<-metaMDS(t(tethys_ordered),autotransform=F)

loping_genera<-unique(loping$matched_name) #creates list of genera found globally in Lopingian

iran_guad_gen<-unique(subset(iran$matched_name,iran$early_interval=="Wordian" | iran$early_interval=="Capitanian")) #creates list of genera found in Guadalupian and Lopingian of each region
iran_lop_gen<-unique(subset(iran$matched_name,iran$early_interval=="Wuchiapingian" | iran$early_interval=="Changhsingian"))
schina_guad_gen<-unique(subset(schina$matched_name,schina$early_interval=="Kungurian" | schina$early_interval=="Roadian" | schina$early_interval=="Wordian" | schina$early_interval=="Capitanian"))
schina_lop_gen<-unique(subset(schina$matched_name,schina$early_interval=="Wuchiapingian" | schina$early_interval=="Changhsingian"))
pakistan_guad_gen<-unique(subset(pakistan$matched_name,pakistan$early_interval=="Wordian" | pakistan$early_interval=="Capitanian"))
pakistan_lop_gen<-unique(subset(pakistan$matched_name,pakistan$early_interval=="Wuchiapingian" | pakistan$early_interval=="Changhsingian"))
tethys_guad_gen<-unique(subset(tethys$matched_name,tethys$early_interval=="Roadian" | tethys$early_interval=="Wordian" | tethys$early_interval=="Capitanian"))
tethys_lop_gen<-unique(subset(tethys$matched_name,tethys$early_interval=="Wuchiapingian" | tethys$early_interval=="Changhsingian"))


#finds "species" scores and identifies as Guadalupian (extinct or extirpated), Lopingian, or both
iran_species<-iran_MDS$species[,c(1:2)]
iran_species<-data.frame(iran_MDS$species[,c(1:2)],"extinct"=ifelse(rownames(iran_species) %in% iran_guad_gen & rownames(iran_species) %in% iran_lop_gen,"Both",ifelse(rownames(iran_species) %in% iran_lop_gen,"Lop",ifelse(rownames(iran_species) %in% loping_genera,"GuadExtirp","GuadExtinct"))))
schina_species<-schina_MDS$species[,c(1:2)]
schina_species<-data.frame(schina_MDS$species[,c(1:2)],"extinct"=ifelse(rownames(schina_species) %in% schina_guad_gen & rownames(schina_species) %in% schina_lop_gen,"Both",ifelse(rownames(schina_species) %in% schina_lop_gen,"Lop",ifelse(rownames(schina_species) %in% loping_genera,"GuadExtirp","GuadExtinct"))))
pakistan_species<-pakistan_MDS$species[,c(1:2)]
pakistan_species<-data.frame(pakistan_MDS$species[,c(1:2)],"extinct"=ifelse(rownames(pakistan_species) %in% pakistan_guad_gen & rownames(pakistan_species) %in% pakistan_lop_gen,"Both",ifelse(rownames(pakistan_species) %in% pakistan_lop_gen,"Lop",ifelse(rownames(pakistan_species) %in% loping_genera,"GuadExtirp","GuadExtinct"))))
tethys_species<-tethys_MDS$species[,c(1:2)]
tethys_species<-data.frame(tethys_MDS$species[,c(1:2)],"extinct"=ifelse(rownames(tethys_species) %in% tethys_guad_gen & rownames(tethys_species) %in% tethys_lop_gen,"Both",ifelse(rownames(tethys_species) %in% tethys_lop_gen,"Lop",ifelse(rownames(tethys_species) %in% loping_genera,"GuadExtirp","GuadExtinct"))))

#To plot Tethys MDS figure (FIG 2)
pdf("fig2_tethys_mds.pdf",width=6,height=6)
plot(tethys_MDS,type="n",ylim=c(min(tethys_species$MDS2),max(tethys_species$MDS2)))
points(tethys_species$MDS1,tethys_species$MDS2,pch=c(2,3,3,4)[tethys_species$extinct],col=c("black","gray","black","black")[tethys_species$extinct])
points(scores(tethys_MDS),pch=c(16,16,16,15,15),cex=c(2.25,2.25,2.25,2,2))
text(scores(tethys_MDS)[,1],scores(tethys_MDS)[,2],strtrim(rownames(scores(tethys_MDS)),1),col="white",cex=0.75)
legend(0.26,0.65,c("Guadalupian","Lopingian"),pch=c(16,15),cex=0.7)
legend(0.56,0.65,c("Lopingian","Both","Guadalupian","Guad Extinct"),pch=c(4,2,3,3),col=c("black","black","black","gray"),cex=0.7)
dev.off()


#To plot regional MDS figures (FIG 4)
pdf("fig4_region_mds.pdf",width=6,height=10.8)
par(mfrow=c(3,1))
par(pin=c(4,2.6))

plot(iran_MDS,type="n")
points(iran_species$MDS1,iran_species$MDS2,pch=c(2,3,3,4)[iran_species$extinct],col=c("black","gray","black","black")[iran_species$extinct])
points(scores(iran_MDS),pch=c(16,16,16,15,15),cex=c(2.25,2.25,2.25,2,2))
text(scores(iran_MDS)[,1],scores(iran_MDS)[,2],strtrim(rownames(scores(iran_MDS)),1),col="white",cex=0.75)
mtext("A. Iran",3,0.3,adj=0)

plot(schina_MDS,type="n")
points(schina_species$MDS1,schina_species$MDS2,pch=c(2,3,3,4)[schina_species$extinct],col=c("black","gray","black","black")[schina_species$extinct])
points(scores(schina_MDS),pch=c(16,16,16,15,15),cex=c(2.25,2.25,2.25,2,2))
text(scores(schina_MDS)[,1],scores(schina_MDS)[,2],strtrim(rownames(scores(schina_MDS)),1),col="white",cex=0.75)
mtext("B. South China",3,0.3,adj=0)

plot(pakistan_MDS,type="n")
points(pakistan_species$MDS1,pakistan_species$MDS2,pch=c(2,3,3,4)[pakistan_species$extinct],col=c("black","gray","black","black")[pakistan_species$extinct])
points(scores(pakistan_MDS),pch=c(16,16,15,15),cex=c(2.25,2.25,2,2))
text(scores(pakistan_MDS)[,1],scores(pakistan_MDS)[,2],strtrim(rownames(scores(pakistan_MDS)),1),col="white",cex=0.75)
legend(-0.73,-0.25,c("Guadalupian","Lopingian"),pch=c(16,15),cex=0.8)
legend(0.44,-0.25,c("Lopingian","Both","Guadalupian","Guad Extinct"),pch=c(4,2,3,3),col=c("black","black","black","gray"),cex=0.75)
mtext("C. Pakistan",3,0.3,adj=0)

dev.off()

