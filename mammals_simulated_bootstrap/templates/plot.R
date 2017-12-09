#!/bin/env Rscript
library(ggplot2)
library(plyr)
library(scales)
library(zoo)
library(plyr)

## FBP
fbp=read.table("!{fbp}",header=T,na.strings = "N/A")
colnames(fbp)=c("tree","brid","length","boot","depth","simu","tbetrue","fbptrue")
fbp=fbp[fbp$depth>1,]
svg("fbp_simu_vs_boot.svg")
	ggplot(fbp,aes(x=simu,y=boot))+geom_point(size=2,alpha=0.2)+ggtitle("FBP")
dev.off()

## TBE
tbe=read.table("!{tbe}",header=T,na.strings = "N/A")
colnames(tbe)=c("tree","brid","length","boot","depth","simu","tbetrue")
tbe=tbe[tbe$depth>1,]

svg("tbe_simu_vs_boot.svg")
	ggplot(tbe,aes(x=simu,y=boot))+geom_point(size=2,alpha=0.2)+ggtitle("TBE")
dev.off()

sink("!{div}_!{seed}_correlations.txt")
print("FBP boot vs. simu")
cor(fbp$boot,fbp$simu)
cor(fbp$simu,fbp$boot,method="spearman")

print("TBE boot vs. simu")
cor(tbe$boot,tbe$simu)
cor(tbe$simu,tbe$boot,method="spearman")

## Comparing to true trees
svg("tbe_true_vs_boot.svg")
ggplot(tbe,aes(x=tbetrue,y=boot))+geom_point(size=2,alpha=0.2)+ggtitle("TBE")
##plot(tbe$tbetrue,tbe$boot,cex=0.5,pch=20, col=alpha("black",0.2))
dev.off()
print("TBE true vs. TBE boot")
cor(tbe$tbetrue,tbe$boot,use="complete.obs")
cor(tbe$tbetrue,tbe$boot,method="spearman")

svg("tbe_true_vs_simu.svg")
plot(tbe$tbetrue,tbe$simu,cex=0.5,pch=20)
dev.off()
print("TBE true vs. TBE simu")
cor(tbe$tbetrue,tbe$simu,use="complete.obs")
cor(tbe$tbetrue,tbe$simu,method="spearman")

svg("tbe_true_vs_fbp_boot.svg")
plot(fbp$tbetrue,fbp$boot,cex=0.5,pch=20)
dev.off()
print("TBE true vs.  FBP boot")
cor(fbp$tbetrue,fbp$boot,use="complete.obs")
cor(fbp$tbetrue,fbp$tbetrue,method="spearman")

svg("tbe_true_vs_fbp_simu.svg")
plot(fbp$tbetrue,fbp$simu,cex=0.5,pch=20)
dev.off()
print("TBE true vs.  FBP simu")
cor(fbp$tbetrue,fbp$simu,use="complete.obs")
cor(fbp$tbetrue,fbp$simu,method="spearman")

print("FBP True vs. FBP Simu")
cor(fbp$fbptrue,fbp$simu,use="complete.obs")
cor(fbp$fbptrue,fbp$simu,method="spearman")

print("FBP True vs. FBP Boot")
cor(fbp$fbptrue,fbp$boot,use="complete.obs")
cor(fbp$fbptrue,fbp$boot,method="spearman")
sink()

svg("fbp_true_vs_fbp_simu.svg")
ggplot(fbp,aes(x=fbptrue,y=simu,color=as.factor(fbptrue)))+geom_boxplot()
dev.off()
svg("fbp_true_vs_fbp_boot.svg")
ggplot(fbp,aes(x=fbptrue,y=boot,color=as.factor(fbptrue)))+geom_boxplot()
dev.off()

## Rogue analysis
tryCatch(
{
        rogues=read.table(gzfile("!{refrogues}"),stringsAsFactors=F)
        rogueboots=read.table(text=paste0(head(readLines("!{tbelogboot}"), -1)),stringsAsFactors=F,sep=":",skip = 8,header=F)
        colnames(rogueboots)=c("name","index")
        roguesimu=read.table(text=paste0(head(readLines("!{tbelogsimu}"), -1)),stringsAsFactors=F,sep=":",skip = 8,header=F)
        colnames(roguesimu)=c("name","index")
        roguesimu$name=gsub(" ","",roguesimu$name)
        rogueboots$name=gsub(" ","",rogueboots$name)
        
        svg("!{div}_!{seed}_rogues_simu_vs_boot.svg")
        	plot(roguesimu$index/10,rogueboots$index/10,
        		pch=20,
        		col=ifelse(rogueboots$name%in%rogues$V1,"orange","blue"),
        		cex=ifelse(rogueboots$name%in%rogues$V1,2,0.5),
        		xlab="Instablity score / Simulated",
        		ylab="Instability score / Bootstrap")
        	legend("topleft", legend=c("Rogues", "Others"),col=c("orange", "blue"), pch=c(20,20),pt.cex=c(2,0.5))
        	abline(v=arrange(roguesimu,desc(index))[100,"index"]/10)
        	abline(h=arrange(rogueboots,desc(index))[100,"index"]/10)
        dev.off()

	# "roc" curves
	nrogues=length(rogues$V1)
	totaltaxa=length(rogueboots$index)
	simucumu=unlist(lapply(1:totaltaxa,function(c){sum(head(arrange(roguesimu,desc(index)), n = c)[,"name"]%in%rogues$V1)}))/nrogues
	bootcumu=unlist(lapply(1:totaltaxa,function(c){sum(head(arrange(rogueboots,desc(index)), n = c)[,"name"]%in%rogues$V1)}))/nrogues
	svg("rogues_cumu.svg",width=7,height=7.5)
	plot(1:totaltaxa,simucumu,type='l',ylab ="% of total rogues found in first x instable taxa",xlab="First x instable taxa",lwd=2)
	lines(1:totaltaxa,bootcumu,col="blue",lwd=2)
	abline(v=100)
	legend("bottomright", legend=c("Bootstrap trees", "Simulation-based trees"),col=c("blue", "black"),lty=1,lwd=2)
	dev.off()
       
        sink("!{div}_!{seed}_correlations.txt",append=T)
		print("Total Taxa")
		print(totaltaxa)
		print("Nb Rogues")
		print(nrogues)
        	print("Number of rogues in the first 100 instable taxa (simu)")
		print(table(head(arrange(roguesimu,desc(index)), n = 100)[,"name"]%in%rogues$V1))
        	print("Number of rogues in the first 100 instable taxa (boot)")
        	print(table(head(arrange(rogueboots,desc(index)), n = 100)[,"name"]%in%rogues$V1))
        sink()
   },
   error=function(cond) {
            return(NA)
   },
   warning=function(cond) {
            return(NULL)
   },
   finally={
   }
)
