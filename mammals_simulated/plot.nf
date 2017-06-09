params.tbetree="results/indelible/supports/1_65512081_tbe.nw"
params.fbptree="results/indelible/supports/1_65512081_fbp.nw"
params.noisytbetree="results/indelible0.5/supports/1_65512081_tbe.nw"
params.noisyfbptree="results/indelible0.5/supports/1_65512081_fbp.nw"
params.resultdir="results/plots"
params.allsupports="results/indelible/allsupports.txt"
params.conflicts="results/indelible/quartets/groupedinfos.txt"
params.noisyallsupports="results/indelible0.5/allsupports.txt"
params.noisyconflicts="results/indelible0.5/quartets/groupedinfos.txt"

params.boosterlog="results/indelible0.5/supports/1_65512081_tbe.log"
params.roguefile="results/indelible0.5/aligns/rogues.txt"

tbetree=file(params.tbetree)
fbptree=file(params.fbptree)
noisytbetree=file(params.noisytbetree)
noisyfbptree=file(params.noisyfbptree)

resultdir=file(params.resultdir)
resultdir.with{mkdirs()}

rawallsupports=file(params.allsupports)
rawconflicts=file(params.conflicts)
noisyallsupports=file(params.noisyallsupports)
noisyconflicts=file(params.noisyconflicts)

boosterlog=file(params.boosterlog)
roguefile=file(params.roguefile)

// Depth Histograms for non noisy trees
process FigureS10Depth{
	
	input:
	file allsupport from rawallsupports
	each cutoff from 0.7

	output:
	file "*.svg" into figs10plots1

	shell:
	outprefix="FigureS10_no_noise"
	template 'plot_depth_histograms.R'
}

// Size Histograms for non noisy tree
process FigureS10Size{

	input:
	file allsupport from rawallsupports
	each cutoff from 0.7

	output:
	file "*.svg" into figs10plots2

	shell:
	outprefix="FigureS10_no_noise"
	template 'plot_size_histograms.R'
}

// Conflict Histograms for non noisy tree
process FigureS10Conflict{

	input:
	file conflict from rawconflicts
	each cutoff from 0.7

	output:
	file "*.svg" into figs10plots3

	shell:
	outprefix="FigureS10_no_noise"
	template 'plot_conflict_histograms.R'
}


figs10plots1.mix(figs10plots2,figs10plots3).subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

// Depth Histograms for noisy tree
process FigureS10S11Depth{

	input:
	file allsupport from noisyallsupports
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs10s11plots1

	shell:
	outprefix="FigureS10S11_noisy"
	template 'plot_depth_histograms.R'
}

// Size Histograms for RAxML tree
process FigureS10S11Size{

	input:
	file allsupport from noisyallsupports
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs10s11plots2

	shell:
	outprefix="FigureS10S11_noisy"
	template 'plot_size_histograms.R'
}

// Conflict Histograms for RAxML tree
process FigureS10S11Conflict{

	input:
	file conflict from noisyconflicts
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs10s11plots3

	shell:
	outprefix="FigureS10S11_noisy"
	template 'plot_conflict_histograms.R'
}

figs10s11plots1.mix(figs10s11plots2,figs10s11plots3).subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

process FigureS12Rogues{

	input:
	file boosterlog
	file roguefile

	output:
	file "*.svg" into rogueplots mode flatten

	shell:
	'''
	#!/usr/bin/env Rscript
	library(ggplot2)
	cbPalette <- c("#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")
	moving=read.table("!{boosterlog}",skip=8,nrow=length(readLines("!{boosterlog}")) - 9, header=F)
	moving=moving[,c(1,3)]
	colnames(moving)=c("taxon","tindex")
	rogue=read.table("!{roguefile}",header=F)
	moving$rogue=ifelse(moving$taxon%in%rogue$V1,"Rogues","Others")

	svg("FigureS12_rogue_transfer_index.svg",width=8,height=7)
	print(ggplot(moving,aes(x=factor(rogue),tindex))+
	  geom_boxplot()+
	  theme_bw()+
	  scale_fill_manual(values=cbPalette)+
	  scale_colour_manual(values=cbPalette)+
	  ylab("Instability Index")+
	  xlab(""))
	dev.off()
	'''
}

rogueplots.subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}
