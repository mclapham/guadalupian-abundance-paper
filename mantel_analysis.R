#MUST HAVE VEGAN PACKAGE DOWNLOADED
#RUN OCCURRENCE_ACQ.R FIRST
#calculates dissimilarity matrix of stage pairs based on faunal composition
#compares to dissimilarity matrix of stage pairs based on the geographic distribution of sampling
#uses mantel test to assess correlation of dissimilarity matrices

library(vegan)

#calculates Bray-Curtis distance between each stage pair based on genus occurrence frequencies
stage_dist<-vegdist(t(tethys_ordered))

mantel_results<-data.frame(bin=numeric(0),stat=numeric(0),p.val=numeric(0))

for (n in seq(0.5,5,by=0.5)) {
  #groups occurrences into n by n degree bins
  occs<-tethys
  
  #rounds paleolatitude
  occs$paleolat<-round(occs$paleolat/n)*n
  occs$paleolng<-round(occs$paleolng/n)*n
  
  #assigns bin number to each occurrence
  occs$bin_no<-as.factor(paste(occs$paleolng,occs$paleolat))
  
  #applies occurrence counting function to each stage
  stage_occs<-sapply(split(occs,occs$cx_int_no),function(x) sapply(split(x$matched_name,x$bin_no),length))
  
  #puts in chronological order
  stage_occs<-stage_occs[,order(match(as.numeric(colnames(stage_occs)),time_int$interval_no),decreasing=T)]
  
  #renames columns with interval name
  colnames(stage_occs)<-time_int$interval_name[match(as.numeric(colnames(stage_occs)),time_int$interval_no)]
  
  #converts counts to presence/absence data
  stage_occs_pres<-apply(stage_occs,2,function(x) ceiling(x/sum(x)))
  
  #calculates Jaccard distance between all stage pairs based on presence/absence of samples in lat/long bins
  geog_occs_dist<-vegdist(t(stage_occs_pres),method="jaccard")
  
  #performs Mantel test to calculate correlation between sampling-distance matrix and taxonomic-distance matrix
  mantel_temp<-mantel(geog_occs_dist,stage_dist,method="spearman",permutations=5000)
  mantel_results<-rbind(mantel_results,data.frame(bin=n,stat=mantel_temp$statistic,p.val=mantel_temp$signif))
}

pdf("fig3_mantel_results.pdf",width=6,height=6)
par(mar=c(5,4,4,4))
plot(mantel_results$bin,mantel_results$stat,ylim=c(0,1),xlab="Geographic bin size (degrees)",ylab=expression(paste("Correlation coefficient (",rho,")")),type="o",lwd=1.5,pch=16,cex=1.5)
lines(mantel_results$bin,mantel_results$p.val*2,lty=2,lwd=1.5)
points(mantel_results$bin,mantel_results$p.val*2,pch=22,bg="white",cex=1.5)
axis(side=4,at=seq(0,0.5,by=0.05)*2,labels=seq(0,0.5,by=0.05))
abline(h=0.1,lty=3)
mtext("P value",side=4,line=3)
legend(0.5,1,c("Correlation","P value"),cex=0.9,lty=c(1,2),pch=c(16,22),pt.bg="white",bty="n")
dev.off()
