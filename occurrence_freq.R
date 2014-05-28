#CODE TO DOWNLOAD BRACHIOPOD, BIVALVE, GASTROPOD OCCURRENCES
#CALCULATES OCCURRENCE FREQUENCY (PROPORTION OF OCCURRENCES) FOR THE THREE GROUPS
#PLOTS OCCURRENCE FREQUENCY PER STAGE WITH BINOMIAL CONFIDENCE INTERVALS

#Specifies input parameters
include_taxon<-"Rhynchonelliformea, Bivalvia, Gastropoda"
maxinterval<-"Norian"
mininterval<-"Kimmeridgian"

#reads names and ages of time intervals
time_int<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=all&limit=all")

#finds maximum age of oldest user-specified interval
max_interval_ma<-subset(time_int$early_age,time_int$interval_name==maxinterval)

#finds minimum age of youngest user-specified interval
min_interval_ma<-subset(time_int$late_age,time_int$interval_name==mininterval)

#reads occurrences based on specified taxa and age range
bbg_occurrences<-read.csv(paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",include_taxon,"&min_ma=",min_interval_ma,"&max_ma=",max_interval_ma,"&show=time,geo,phylo,ident&limit=all",sep=""))

#finds only those collections resolved to a stage
resolved_bbgoccs<-subset(bbg_occurrences,bbg_occurrences$cx_int_no %in% subset(time_int$interval_no,time_int$level==5))

#deletes occurrences not resolved to at least genus level using classified and unclassified species
resolved_bbgoccs<-subset(resolved_bbgoccs,resolved_bbgoccs$matched_rank<=5)

#deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
resolved_bbgoccs<-subset(resolved_bbgoccs,resolved_bbgoccs$genus_reso=="" | resolved_bbgoccs$genus_reso=="n. gen.")

#list of marine environments
carbonate_env<-c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")
siliciclastic_env<-c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
marine_env<-c(carbonate_env,siliciclastic_env)

#finds collections with appropriate environments
resolved_bbgoccs<-subset(resolved_bbgoccs,resolved_bbgoccs$environment %in% marine_env)

#extracts just the genus name from character string
resolved_bbgoccs$matched_name<-gsub(" .*$","",resolved_bbgoccs$matched_name)

#function to calculate binomial confidence intervals, given probability p, sample size n, and confidence level alpha
binconf<-function(p,n,alpha) {
  # Estimate based on Wilson score interval
  alpha1<-(alpha+1)/2
  z<-qnorm(alpha1)
  maxconf<<-(p+(z^2)/(2*n) + z * sqrt(p*(1-p)/n+(z^2)/(4*n^2)))/(1+(z^2)/n)
  minconf<<-(p+(z^2)/(2*n) - z * sqrt(p*(1-p)/n+(z^2)/(4*n^2)))/(1+(z^2)/n)
  conf<<-maxconf-p
  return(list(maxconf=maxconf,minconf=minconf))
  }

#calculates number of occurrences per stage
brach_occs<-sapply(split(resolved_bbgoccs,resolved_bbgoccs$cx_int_no),function(x) nrow(subset(x,x$class=="Strophomenata" | x$class=="Rhynchonellata" | x$class=="Chileata")))
biv_occs<-sapply(split(resolved_bbgoccs,resolved_bbgoccs$cx_int_no),function(x) nrow(subset(x,x$class=="Bivalvia")))
gast_occs<-sapply(split(resolved_bbgoccs,resolved_bbgoccs$cx_int_no),function(x) nrow(subset(x,x$class=="Gastropoda")))
all_occs<-sapply(split(resolved_bbgoccs,resolved_bbgoccs$cx_int_no),nrow)

#calculates occurrence frequency
brach_prop<-brach_occs/all_occs
biv_prop<-biv_occs/all_occs
gast_prop<-gast_occs/all_occs

#puts in chronological order
brach_prop<-brach_prop[order(match(as.numeric(names(brach_prop)),time_int$interval_no),decreasing=T)]
biv_prop<-biv_prop[order(match(as.numeric(names(biv_prop)),time_int$interval_no),decreasing=T)]
gast_prop<-gast_prop[order(match(as.numeric(names(gast_prop)),time_int$interval_no),decreasing=T)]
all_occs<-all_occs[order(match(as.numeric(names(all_occs)),time_int$interval_no),decreasing=T)]

#renames columns with interval name
names(brach_prop)<-time_int$interval_name[match(as.numeric(names(brach_prop)),time_int$interval_no)]
names(biv_prop)<-time_int$interval_name[match(as.numeric(names(biv_prop)),time_int$interval_no)]
names(gast_prop)<-time_int$interval_name[match(as.numeric(names(gast_prop)),time_int$interval_no)]

#calculates binomial confidence intervals
brach_conf<-binconf(brach_prop,all_occs,0.95)
biv_conf<-binconf(biv_prop,all_occs,0.95)
gast_conf<-binconf(gast_prop,all_occs,0.95)

#plots results
pdf("brach_gast_occs.pdf",width=10.5)
plot(brach_prop,pch=16,cex=1.5,ylim=c(0,1),xlab="",ylab="Occurrence frequency",xaxt="n")

#adds 95% binomial confidence intervals for proportions
segments(seq(1,length(brach_prop)),brach_conf$maxconf,seq(1,length(brach_prop)),brach_conf$minconf,lwd=1.5)
segments(seq(1,length(biv_prop)),biv_conf$maxconf,seq(1,length(biv_prop)),biv_conf$minconf,lwd=1.5)
segments(seq(1,length(gast_prop)),gast_conf$maxconf,seq(1,length(gast_prop)),gast_conf$minconf,lwd=1.5)

#adds points for bivalve and gastropod occurrence frequency
points(biv_prop,cex=1.5,pch=21,bg="gray")
points(gast_prop,cex=1.5,pch=21,bg="white")

par(xpd=T)
rect(seq(0.52,12.48,length.out=13),-0.04,(seq(0.52,12.48,length.out=13)+0.9975),-0.1)
par(xpd=F)
mtext(c("Ass","Sakm","Art","Kung","Road","Word","Cap","Wuch","Chang","Ind","Olen","Anis","Lad"),1,0.15,at=seq(1,13))
mtext("Stage",1,2.5)

#add lines for Guadalupian and end-Permian extinctions
abline(v=7.5,lty=2)
abline(v=9.49)

legend(10.2,1.05,c("Rhynchonelliformea","Bivalvia","Gastropoda"),pch=c(16,21,21),pt.bg=c("black","gray","white"),cex=1.2,bty="n")
dev.off()
#end of plot
