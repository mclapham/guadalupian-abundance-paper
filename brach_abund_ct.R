#read data file with census counts
abund<-read.csv("https://github.com/mclapham/guadalupian-abundance-paper/blob/master/pt_abund_data.csv?raw=T")

#changes name from Early/Lower Triassic to Early Triassic
levels(abund$epoch)<-c("Cisuralian","Early Triassic","Guadalupian","Lopingian")

#function to calculate prportional abundance of rhynchonelliform brachiopods per collection
brach.ct<-function(coll_data) {
  sapply(split(coll_data,coll_data$collection_no),function(x) sum(subset(x$abund_value,x$class_name=="Strophomenata" | x$class_name=="Rhynchonellata" | x$class_name=="Chileata"))/sum(x$abund_value))
}

#calculates proportional abundance per epoch
brach_abund<-sapply(split(abund,abund$epoch),brach.ct)

pdf("brach_abund.pdf",width=7)
boxplot(brach_abund[c(1,3,4,2)],ylab="Rhynchonelliformea proportional abundance")
dev.off()
