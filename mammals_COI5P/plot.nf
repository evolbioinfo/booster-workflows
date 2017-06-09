params.tbetreefasttree="results/fasttree/support/1_31144_tbe.nw"
params.fbptreefasttree="results/fasttree/support/1_31144_fbp.nw"
params.namemapfasttree="results/fasttree/aligns/name_map.txt.gz"

params.tbetreeraxml="results/raxml/support/1_31144_tbe.nw"
params.fbptreeraxml="results/raxml/support/1_31144_fbp.nw"
params.namemapraxml="results/raxml/aligns/name_map.txt.gz"

params.resultdir="results/plots"
params.ncbi="../data/ncbitax/ncbi_labels.nw"

params.allsupportfasttree="results/fasttree/allsupports.txt"
params.allsupportraxml="results/raxml/allsupports.txt"
params.conflictfasttree="results/fasttree/quartets/groupedinfos.txt"
params.conflictraxml="results/raxml/quartets/groupedinfos.txt"

params.raxmltree="results/raxml/trees/ref_1_31144.nw.gz"
params.raxmlboottrees="results/raxml/trees/boot_1_31144.nw.gz"
params.fasttreetree="results/fasttree/trees/ref_1_31144.nw.gz"
params.fasttreeboottrees="results/fasttree/trees/boot_1_31144.nw.gz"

tbetreefasttree=file(params.tbetreefasttree)
fbptreefasttree=file(params.fbptreefasttree)
namemapfasttree=file(params.namemapfasttree)
tbetreeraxml=file(params.tbetreeraxml)
fbptreeraxml=file(params.fbptreeraxml)
namemapraxml=file(params.namemapraxml)
ncbi=file(params.ncbi)
resultdir=file(params.resultdir)
resultdir.with{mkdirs()}

allsupportfasttree=file(params.allsupportfasttree)
allsupportraxml=file(params.allsupportraxml)
conflictfasttree=file(params.conflictfasttree)
conflictraxml=file(params.conflictraxml)

raxmltree=file(params.raxmltree)
raxmlboottrees=file(params.raxmlboottrees)
fasttreetree=file(params.fasttreetree)
fasttreeboottrees=file(params.fasttreeboottrees)

// We plot Figure 2 Right Panel
// The left panel will be ploted
// In FigS3S4 and simulated data
process Figure2Right{
	input:
	file fbptree from fbptreefasttree
	file tbetree from tbetreefasttree
	file namemap from namemapfasttree
	file ncbi

	output:
	file "*.svg" into fig2plots mode flatten

	shell:
	outprefix="Figure2_Right"
	template 'draw_simians.sh'
}

fig2plots.subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

// We plot Figure S5 (simian tree with RAxML)
process FigureS5{

	input:
	file fbptree from fbptreeraxml
	file tbetree from tbetreeraxml
	file namemap from namemapraxml
	file ncbi

	output:
	file "*.svg" into figs5plots mode flatten

	shell:
	outprefix="FigureS5"
	template 'draw_simians.sh'
}

figs5plots.subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

// We count the number of times Maxomys rajah and Canis adustus
// must be moved around real simian clade or around inferred simian
// clade (1 error)
process FigS5RAxMLStats {

	cpus 12

	input:
	file tree from raxmltree
	file boottrees from raxmlboottrees
	file namemap from namemapraxml
	file ncbi

	output:
	stdout into figs5statsraxml

	shell:
	template 'taxa_transfer_stats.sh'
}

figs5statsraxml.collectFile(name: 'figs5statsraxml.txt').subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

// We count the number of times Maxomys rajah and Canis adustus
// must be moved around real simian clade or around inferred simian
// clade (2 errors)
process FigS5FastTreeStats {

	cpus 12

	input:
	file tree from fasttreetree
	file boottrees from fasttreeboottrees
	file namemap from namemapfasttree
	file ncbi

	output:
	stdout into figs5statsfasttree

	shell:
	template 'taxa_transfer_stats.sh'
}

figs5statsfasttree.collectFile(name: 'figs5statsfasttree.txt').subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

// Depth Histograms for FastTree tree
process FigureS3Depth{
	
	input:
	file allsupport from allsupportfasttree
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs3plots1

	shell:
	outprefix="FigureS3"
	template 'plot_depth_histograms.R'
}

// Size Histograms for FastTree tree
process FigureS3Size{

	input:
	file allsupport from allsupportfasttree
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs3plots2

	shell:
	outprefix="FigureS3"
	template 'plot_size_histograms.R'
}

// Conflict Histograms for FastTree tree
process FigureS3Conflict{

	input:
	file conflict from conflictfasttree
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs3plots3

	shell:
	outprefix="FigureS3"
	template 'plot_conflict_histograms.R'
}


figs3plots1.mix(figs3plots2,figs3plots3).subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

// Depth Histograms for RAxML tree
process FigureS4Depth{

	input:
	file allsupport from allsupportraxml
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs4plots1

	shell:
	outprefix="FigureS4"
	template 'plot_depth_histograms.R'
}

// Size Histograms for RAxML tree
process FigureS4Size{

	input:
	file allsupport from allsupportraxml
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs4plots2

	shell:
	outprefix="FigureS4"
	template 'plot_size_histograms.R'
}

// Conflict Histograms for RAxML tree
process FigureS4Conflict{

	input:
	file conflict from conflictraxml
	each cutoff from 0.5,0.7,0.9

	output:
	file "*.svg" into figs4plots3

	shell:
	outprefix="FigureS4"
	template 'plot_conflict_histograms.R'
}


figs4plots1.mix(figs4plots2,figs4plots3).subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}

// Figure S6: Full Mammalian RAxML tree with annotations
// and TBE&FBP supported branches
process FigureS6{
	input:
	file fbptree from fbptreeraxml
	file tbetree from tbetreeraxml
	file namemap from namemapraxml
	file ncbi

	output:
	file "*.svg" into figs6plots mode flatten

	shell:
	outprefix="FigureS6"
	template 'draw_full.sh'
}

// Figure S6: Stats about Figure S6
// Getting all NCBI clades and their TBE and FBP support in full
// RAxML trees
process FigureS6Stats{
	input:
	file fbptree from fbptreeraxml
	file tbetree from tbetreeraxml
	file namemap from namemapraxml
	file ncbi

	output:
	file "*_supports.txt" into figs6stats mode flatten

	shell:
	'''
	#!/usr/bin/env bash
	gotree rename -i !{tbetree} -m !{namemap} -r > tmp
	gotree rename -i !{fbptree} -m !{namemap} -r > tmpfbp
	gotree unroot -i !{ncbi} > ncbi_unroot

	gotree prune -i tmp    -c !{ncbi} | gotree compare edges -m --moved-taxa -c ncbi_unroot | awk -F "\t" '{if($7>1){print $0}}' > ncbi_compare
	gotree prune -i tmpfbp -c !{ncbi} | gotree compare edges -m --moved-taxa -c ncbi_unroot | awk -F "\t" '{if($7>1){print $0}}' > ncbi_comparefbp

	tail -n+2 ncbi_compare          | sort -u -k12 -k10,10n -t$'\t' | awk -F"\t" '!_[$12]++' > ncbi_tbe_supports.txt
	tail -n+2 ncbi_comparefbp       | sort -u -k12 -k10,10n -t$'\t' | awk -F"\t" '!_[$12]++' > ncbi_fbp_supports.txt
	'''
}

figs6plots.mix(figs6stats).subscribe{
	f->f.copyTo(resultdir.resolve(f.name))
}
