params.tbetree="results/supports/1_691607951_tbe.nw"
params.fbptree="results/supports/1_691607951_fbp.nw"
params.namemap="results/aligns/name_map.txt.gz"
params.resultdir="results/plots"
params.recomb="results/formated.txt"

params.boosterlog="results/supports/1_691607951_tbe.log"
params.allsupports="results/allsupports.txt"

tbetree=file(params.tbetree)
fbptree=file(params.fbptree)
namemap=file(params.namemap)
resultdir=file(params.resultdir)
resultdir.with{mkdirs()}

allsupports=file(params.allsupports)
recomb=file(params.recomb)
boosterlog=file(params.boosterlog)

mediumtrees=Channel.from(["Sample1", file("results/supportstrees/16_1339632896_fbp.nw"),file("results/supportstrees/16_1339632896_tbe.nw")],
			 ["Sample2", file("results/supportstrees/16_1374534207_fbp.nw"),file("results/supportstrees/16_1374534207_tbe.nw")])

mediumtrees.into{mediumtrees1; mediumtrees2; mediumtrees3}

process WriteItolAnnotationFile {
	input:
	file recomb

	output:
	file "annot_itol.txt" into itolannot

	shell:
	template 'write_itol_annotation_file.R'

}

process WriteItolDrawingConfigFiles {
	input:
	file recomb

	output:
	set file("config_tbe.txt"),file("config_fbp.txt") into itolconfig

	shell:
	template 'write_itol_config_file.sh'
}

process Figure1Trees{
	input:
	file tbetree
	file fbptree
	file itolannot
	file namemap
	set file(configtbe),file(configfbp) from itolconfig

	output:
	file "*.svg" into fig1plots1 mode flatten

	shell:
	outprefix="Figure1"
	template 'draw_full_hiv_trees.sh'
}

process Figure1Histograms{
	input:
	file allsupports

	output:
	file "*.svg" into fig1plots2 mode flatten

	shell:
	template 'draw_depth_histograms.R'
}

process Figure1SubtypeSupports{
	input:
	file tbetree
	file fbptree
	file recomb
	file namemap

	output:
	file "*.txt" into fig1plots3 mode flatten

	shell:
	outprefix="Figure1"
	template 'compute_subtypes_supports.sh'
}

process Figure1SubtypeGeographicalSupports{
	input:
	file tbetree
	file fbptree
	file recomb
	file namemap

	output:
	file "*.txt" into fig1plots4 mode flatten

	shell:
	outprefix="Figure1"
	template 'compute_subtypes_geo_supports.sh'
}

fig1plots1.mix(fig1plots2,fig1plots3,fig1plots4).subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

process FigureS8Trees{
	input:
	set val(sample), file(fbptree), file(tbetree) from mediumtrees1
	file itolannot
	file namemap
	set file(configtbe),file(configfbp) from itolconfig

	output:
	file "*.svg" into figs8plots1 mode flatten

	shell:
	outprefix="FigureS8_"+sample
	template 'draw_full_hiv_trees.sh'
}

process FigureS8SubtypeSupports{
	input:
	set val(sample), file(fbptree), file(tbetree) from mediumtrees2
	file recomb
	file namemap

	output:
	file "*.txt" into figs8plots2 mode flatten

	shell:
	outprefix="FigureS8_"+sample
	template 'compute_subtypes_supports.sh'
}

process FigureS8SubtypeGeographicalSupports{
	input:
	set val(sample), file(fbptree), file(tbetree) from mediumtrees3
	file recomb
	file namemap

	output:
	file "*.txt" into figs8plots3 mode flatten

	shell:
	outprefix="FigureS8_"+sample
	template 'compute_subtypes_geo_supports.sh'
}


figs8plots1.mix(figs8plots2,figs8plots3).subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}


// Depth Histograms
process FigureS7Depth{
	
	input:
	file allsupport from allsupports
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs7plots1

	shell:
	outprefix="FigureS7"
	template 'plot_depth_histograms.R'
}

// Size Histograms for FastTree tree
process FigureS7Size{

	input:
	file allsupport from allsupports
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs7plots2

	shell:
	outprefix="FigureS7"
	template 'plot_size_histograms.R'
}

figs7plots1.mix(figs7plots2).subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

process FigureS9Recomb{

	input:
	file boosterlog
	file recomb
	file namemap

	output:
	file "*.svg" into figs9plots mode flatten

	shell:
	'''
	#!/usr/bin/env Rscript
	library(ggplot2)
	cbPalette <- c("#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73")
	moving=read.table("!{boosterlog}",skip=8,nrow=length(readLines("!{boosterlog}")) - 9, header=F)
	moving=moving[,c(1,3)]
	colnames(moving)=c("taxon","tindex")
	map=read.table(gzfile("!{namemap}"),header=F)
	rogue=read.table("!{recomb}",header=F)
	rogue=rogue[grep(",",rogue[,2]),1]
	rogue=gsub("_",".",rogue)
	rogue=map[map$V1%in%rogue,2]
	moving$rogue=ifelse(moving$taxon%in%rogue,paste0("Recombinant (",sum(moving$taxon%in%rogue),")"),"Non-recombinant")

	svg("FigureS9_recombinant_transfer_index.svg",width=8,height=7)
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

figs9plots.subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}
