#!/bin/env Rscript
library(ggplot2)
library(plyr)

## FBP
fbp=read.table("!{fbp}",header=T,na.strings = "N/A")
colnames(fbp)=c("tree","brid","length","boot","depth","simu")
fbp=fbp[fbp$depth>1,]
svg("fbp_simu_vs_boot.svg")
	ggplot(fbp,aes(x=simu,y=boot))+geom_point(size=2,alpha=0.2)+ggtitle("FBP")
dev.off()

## TBE
tbe=read.table("!{tbe}",header=T,na.strings = "N/A")
colnames(tbe)=c("tree","brid","length","boot","depth","simu")
tbe=tbe[tbe$depth>1,]

svg("tbe_simu_vs_boot.svg")
	ggplot(tbe,aes(x=simu,y=boot))+geom_point(size=2,alpha=0.2)+ggtitle("TBE")
dev.off()

sink("!{div}_!{seed}_correlations.txt")
	print("FBP boot vs. simu Spearman")
	cor(fbp$boot,fbp$simu,method = "spearman")
	print("FBP boot vs. simu Pearson")
	cor(fbp$boot,fbp$simu,method = "pearson")
	print("TBE boot vs. simu Spearman")
	cor(tbe$boot,tbe$simu,method = "spearman")
	print("TBE boot vs. simu Pearson")
	cor(tbe$boot,tbe$simu,method = "pearson")
sink()	

## Rogue analysis
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

sink("!{div}_!{seed}_correlations.txt",append=T)
	print("Number of rogues in the first 100 instable taxa (simu)")
	table(head(arrange(roguesimu,desc(index)), n = 100)[,"name"]%in%rogues$V1)
	print("Number of rogues in the first 100 instable taxa (boot)")
	table(head(arrange(rogueboots,desc(index)), n = 100)[,"name"]%in%rogues$V1)
sink()	
