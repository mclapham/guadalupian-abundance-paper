#SCRIPT TO ACQUIRE AND FILTER BRACHIOPOD OCCURRENCE DATA
#GENERATES MATRIX OF OCCURRENCE FREQUENCY (COLUMNS=TIME INTERVALS, ROWS=GENERA)
#ORDERED BY TIME INTERVAL FOR TETHYS AND EACH REGION (IRAN, S CHINA, PAKISTAN)

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
occurrences<-read.csv(paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",include_taxon,"&min_ma=",min_interval_ma,"&max_ma=",max_interval_ma,"&show=coords,time,loc,paleoloc,ident&limit=all",sep=""))

#finds only those collections resolved to a stage
resolved_occs<-subset(occurrences,occurrences$cx_int_no %in% subset(time_int$interval_no,time_int$level==5))

#deletes occurrences not resolved to at least genus level using classified and unclassified species
resolved_occs<-subset(resolved_occs,resolved_occs$matched_rank<=5)

#deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
resolved_occs<-subset(resolved_occs,resolved_occs$genus_reso=="" | resolved_occs$genus_reso=="n. gen.")

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
  unsorted_data<-unsorted_data[,order(match(as.numeric(colnames(unsorted_data)),time_int$interval_no),decreasing=T)]
  
  #renames columns with interval name
  colnames(unsorted_data)<-time_int$interval_name[match(as.numeric(colnames(unsorted_data)),time_int$interval_no)]
  
  #removes empty rows
  subset(unsorted_data,apply(unsorted_data,1,function(x) sum(x,na.rm=T))>0)
  
}

#Applies function to re-order columns to chronological order
iran_ordered<-interval.ordering(iran_abund)
schina_ordered<-interval.ordering(schina_abund)
pakistan_ordered<-interval.ordering(pakistan_abund)
tethys_ordered<-interval.ordering(tethys_abund)
