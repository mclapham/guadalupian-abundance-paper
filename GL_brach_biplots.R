#RUN OCCURRENCE_ACQ.R FIRST
#CREATES BIPLOTS OF OCCURRENCE FREQUENCY
#CALCULATES KENDALL RANK ORDER CORRELATION FOR EACH PAIR OF SUCCESSIVE INTERVALS
#PRODUCES FIGS X (Tethys data) and X (individual regions)

#list of all interval names (in order to generate multi-panel plots when some regions are missing data)
int_names<-unique(resolved_occs$early_interval)[order(match(unique(resolved_occs$early_interval),time_int$interval_name),decreasing=T)]

#creates biplots of relative frequency in successive time intervals
abund.biplot<-function(abund_data) {
  for(i in 1:(length(unique(resolved_occs$cx_int_no))-1)) {
      abund_total<-data.frame(rep(0,nrow(abund_data)))
      abund_total[,which(!int_names %in% colnames(abund_data))]<-rep(0,nrow(abund_data)) #adds zero values for time periods not containing data
      abund_total[,which(int_names %in% colnames(abund_data))]<-abund_data #adds data to time periods where it exists
      names(abund_total)<-int_names
      abund_temp<-subset(abund_total[,c(i,i+1)],abund_total[,i]>0 | abund_total[,i+1]>0)
      if(nrow(abund_temp)>0) {
        if(sum(abund_temp[1])>0 & sum(abund_temp[2])>0) {
          plot(abund_temp,pch=16)
          cor_val<-round(cor.test(abund_temp[,1],abund_temp[,2],method="k")$estimate,2)
          p_col<-if(round(cor.test(abund_temp[,1],abund_temp[,2],method="k")$p.value,2)>0.05) {"dark gray"} else {"black"}
          mtext(bquote(tau*" = "*.(cor_val)),3,0.3,adj=1,cex=0.75,col=p_col)
        }  else {
          plot(1,1,type="n",xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),xlab=colnames(abund_temp)[1],ylab=colnames(abund_temp)[2])
          text(0.5,0.5,"No Data")
        }
      } else {
          plot(1,1,type="n",xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),xlab=colnames(abund_temp)[1],ylab=colnames(abund_temp)[2])
          text(0.5,0.5,"No Data")
      }
   }
}

#plots entire Tethys data (figure not included in manuscript)
pdf("figX_tethyan_biplot.pdf",width=13,height=3)
par(mar=c(3,3,2,1.5))
par(mfrow=c(1,length(unique(resolved_occs$cx_int_no))-1))
par(mgp=c(2,0.75,0))
abund.biplot(tethys_ordered)
dev.off()


#plots panels for Iran, South China, and Pakistan (FIG 5)
pdf("fig5_region_biplot.pdf",width=9,height=5)
par(mar=c(3,3,2,1.5))
par(mfrow=c(3,length(unique(resolved_occs$cx_int_no))-1))
par(omi=c(0,0,0,0.3))
par(mgp=c(2,0.75,0))

abund.biplot(iran_ordered)
abund.biplot(schina_ordered)
abund.biplot(pakistan_ordered)

mtext("Iran",side=4,cex=1.1,outer=T,at=0.845,line=0.25)
mtext("South China",side=4,cex=1.1,outer=T,line=0.25)
mtext("Pakistan",side=4,cex=1.1,outer=T,at=0.173,line=0.25)
dev.off()
