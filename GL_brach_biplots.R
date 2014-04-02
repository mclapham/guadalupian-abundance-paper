#Specifies input parameters
include_taxon<-"Brachiopoda"
maxinterval<-"Kungurian"
mininterval<-"Changhsingian"

#reads names and ages of time intervals
time_int<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=all&limit=all")

#finds maximum age of oldest user-specified interval
max_interval_ma<-subset(time_int$early_age,time_int$interval_name==maxinterval)

#finds minimum age of youngest user-specified interval
min_interval_ma<-subset(time_int$late_age,time_int$interval_name==mininterval)

#reads occurrences based on specified taxa and age range
occurrences<-read.csv(paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",include_taxon,"&min_ma=",min_interval_ma,"&max_ma=",max_interval_ma,"&show=coords,time,loc,paleoloc&limit=all",sep=""))

#finds only those collections resolved to a stage
resolved_occs<-subset(occurrences,occurrences$cx_int_no %in% subset(time_int$interval_no,time_int$level==5))

#deletes occurrences not resolved to at least genus level using classified and unclassified species
resolved_occs<-subset(resolved_occs,resolved_occs$taxon_rank=="species" | resolved_occs$taxon_rank=="genus")

#deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
if(length(grep("\\? ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("\\? ",resolved_occs$taxon_name),]
if(length(grep("\" ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("\" ",resolved_occs$taxon_name),]
if(length(grep("cf. ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("cf. ",resolved_occs$taxon_name),]
if(length(grep("aff. ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("aff. ",resolved_occs$taxon_name),]

#extracts just the genus name from character string
resolved_occs$matched_name<-gsub(" .*$","",resolved_occs$matched_name)

#Creates data frames for sub-regions
iran<-subset(resolved_occs,resolved_occs$geoplate==512)
schina<-subset(resolved_occs,resolved_occs$geoplate==611)
pakistan<-subset(resolved_occs,resolved_occs$cc=="PK" & resolved_occs$lat<33.5)
tethys<-subset(resolved_occs,resolved_occs$lng>5 & resolved_occs$lng<145 & resolved_occs$lat>0 & resolved_occs$lat<45)

iran$matched_name<-as.factor(iran$matched_name)
schina$matched_name<-as.factor(schina$matched_name)
pakistan$matched_name<-as.factor(pakistan$matched_name)
tethys$matched_name<-as.factor(tethys$matched_name)

#creates matrix of genus occurrence counts per stage
iran_abund<-sapply(split(iran,iran$cx_int_no),function(x) sapply(split(x$matched_name,x$matched_name),length))
schina_abund<-sapply(split(schina,schina$cx_int_no),function(x) sapply(split(x$matched_name,x$matched_name),length))
pakistan_abund<-sapply(split(pakistan,pakistan$cx_int_no),function(x) sapply(split(x$matched_name,x$matched_name),length))
tethys_abund<-sapply(split(tethys,tethys$cx_int_no),function(x) sapply(split(x$matched_name,x$matched_name),length))

#converts to relative proportion of occurrences
iran_abund<-apply(iran_abund,2,function(x) x/sum(x)) 
schina_abund<-apply(schina_abund,2,function(x) x/sum(x))
pakistan_abund<-apply(pakistan_abund,2,function(x) x/sum(x))
tethys_abund<-apply(tethys_abund,2,function(x) x/sum(x))

#Subfunction to put columns in ascending chronological order
interval.ordering<-function(unsorted_data) {
    #the intervals are in ascending numerical order by interval_no, but not necessarily chronological (this finds the chronological order)
    interval_order<-order(sapply(colnames(unsorted_data),function(x) which(time_int$interval_no==x)),decreasing=T)
                     
    #moves columns into chronological order
    unsorted_data<-unsorted_data[,interval_order]
                     
    #renames columns with interval name
    colnames(unsorted_data)<-time_int$interval_name[sapply(colnames(unsorted_data),function(x) which(time_int$interval_no==x)]
                     
    #removes empty rows
    subset(unsorted_data,apply(unsorted_data,1,function(x) sum(x,na.rm=T))>0)
                     
 }
                   
#Applies function to re-order columns to chronological order
iran_ordered<-interval.ordering(iran_abund)
schina_ordered<-interval.ordering(schina_abund)
pakistan_ordered<-interval.ordering(pakistan_abund)
tethys_ordered<-interval.ordering(tethys_abund)

int_names<-subset(time_int$interval_name,time_int$interval_name %in% unique(resolved_occs$early_interval))[seq(length(unique(resolved_occs$early_interval)),1)]

#creates biplots of relative frequency in successive time intervals

abund.biplot<-function(abund_data) {
  for(i in 1:(length(unique(resolved_occs$cx_int_no))-1)) {
      abund_total<-data.frame(rep(0,nrow(abund_data)))
      abund_total[,which(!int_names %in% colnames(abund_data))]<-rep(0,nrow(abund_data))
      abund_total[,which(int_names %in% colnames(abund_data))]<-abund_data
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


#plots panels for Iran, South China, and Pakistan
pdf("brach_abund_biplot.pdf",width=9,height=5)
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

#plots panels for entire Tethys data
pdf("brach_tethyan_biplot.pdf",width=13,height=3)
par(mar=c(3,3,2,1.5))
par(mfrow=c(1,5))
par(mgp=c(2,0.75,0))
abund.biplot(tethys_ordered)
dev.off()
