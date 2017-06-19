params.size=1449
params.nbref=100
params.nbrand=1000
params.resdir="result"
params.reftree="yuletree"
params.randtree="yuletree"
params.divisions = 10 /* Number of processes per ref tree */

rand = new Random()
refseed   = rand.nextInt()

size = params.size
/* Convert size in tree depth (for balanced trees) */
depth   = (int)(Math.log(size) / Math.log(2))
nbref   = params.nbref
nbrand  = params.nbrand
reftree = params.reftree
randtree= params.randtree

resdir    = file(params.resdir)
resdir.with {mkdirs()}

divisions = params.divisions
nbrandperprocess = (int)(nbrand / params.divisions)

process reftree {

	input:
	val size
	val nbref
	val reftree
	val depth

	output:
	file "reftree*" into reftreechan mode flatten

	shell:
	if( params.reftree == 'balancedtree' )
	'''
	#!/bin/bash
	gotree generate !{reftree} -d !{depth} -n !{nbref} -s $RANDOM \
	       | gotree shuffletips -s $RANDOM \
	       | gotree divide -o reftree
	gzip reftree*
	'''
	else
	'''
	#!/bin/bash
	gotree generate !{reftree} -l !{size} -n !{nbref} -s $RANDOM \
	       | gotree shuffletips -s $RANDOM \
	       | gotree divide -o reftree
	gzip reftree*
	'''
}

process gentrees {

	input:
	val size
	val depth
	val randtree
	val nbrandperprocess
	val divisions
	file reftree from reftreechan

	output:
	set file(reftree), file("randtrees.*.nw.gz") into refrandtreeschan mode flatten

	shell:
	if( params.randtree == 'balancedtree' )
	'''
	#!/bin/bash
	for i in {1..!{divisions}}
	do
		gotree generate !{randtree} -d !{depth} -n !{nbrandperprocess} -s $RANDOM \
	       | gotree shuffletips -s $RANDOM \
	       | gzip -c - \
	       > randtrees.$i.nw.gz
	done
	'''
	else
	'''
	#!/bin/bash
	for i in {1..!{divisions}}
	do
	gotree generate !{randtree} -l !{size} -n !{nbrandperprocess} -s $RANDOM \
	       | gotree shuffletips -s $RANDOM \
	       | gzip -c - \
	       > randtrees.$i.nw.gz
	done
	'''
}

refrand=refrandtreeschan.map{ref,rand -> [rand.name.tokenize('.')[1],ref,rand]}

process comptransfer {

input:
	set val(div), file(ref), file(rand) from refrand

	output:
	set val("${ref.name}"),file("stats_*.txt") into stats

	shell:
	'''
	#!/bin/bash
	gotree compare edges -i !{ref} -c !{rand} -m \
	            | awk '{FS="\t"}{if(NR>1){print ((!{div}-1)*!{nbrandperprocess}+$1) "\t" $2 "\t" $7 "\t" $10}}' \
		    > stats_!{div}_!{params.size}_ref!{params.reftree}_rand!{params.randtree}.txt
	'''
}

groupedstats = stats.groupTuple(by:0)

process ConcatStats {
	input:
	set val(ref), file(statfiles) from groupedstats

	output:
	file "stats_${params.size}_ref${params.reftree}_rand${params.randtree}.txt" into concatstats

	shell:
	'''
	#!/usr/bin/env bash
	cat !{statfiles} > stats_!{params.size}_ref!{params.reftree}_rand!{params.randtree}.txt
	'''
}

process meanDistance{
	input:
	file concatstats

	output:
	file "meandist_*.txt" into meandist

	shell:
	'''
	#!/usr/bin/env Rscript
	read.table("!{concatstats}")->st
	colnames(st)=c("tree","branch","depth","dist")
	aggdata <-aggregate(st$dist, by=list(st$branch,st$depth), FUN=mean)
	write.table(aggdata,"meandist_!{params.size}_ref!{params.reftree}_rand!{params.randtree}.txt",col.names=F, row.names=F, quote=F)
	'''
}

meandist.collectFile(name: 'stats_'+params.size+'_ref'+params.reftree+'_rand'+params.randtree+'.txt').into{groupedstats; groupedstatscopy}

groupedstatscopy.subscribe{
	file -> file.copyTo(resdir.resolve(file.name))
}

process plotstats {

	input:
	file(stats) from groupedstats

	output:
	file("*.svg") into plots mode flatten

	shell:
	'''
	#!/usr/bin/env Rscript
	library(ggplot2)

	q95 <- function(x) {quantile(x,probs=0.95)}

	s=read.table("!{stats}")
	colnames(s)=c("branch","depth","dist")
	svg("stats_!{params.size}_ref!{params.reftree}_rand!{params.randtree}.svg",width=7,height=7)
	ggplot(s, aes(x=depth, y=dist)) +
		  stat_summary(fun.y="mean", geom="point")+
		  ylim(0,max(s$depth))+
		  geom_abline(intercept = -1, slope = 1,color="blue")+
		  theme(text = element_text(size=20))+
		  xlab("Branch min subtree size")+
		  ylab("Transfer distance")+
		  ggtitle("!{size} Taxa")
	dev.off()

	svg("stats_!{params.size}_ref!{params.reftree}_rand!{params.randtree}_norm.svg",width=7,height=3)
	s=s[s$depth>1,]
	s$norm=1-s$dist/(s$depth-1)
	ggplot(s, aes(x=depth, y=norm)) +
		  stat_summary(fun.y=mean,  geom = 'ribbon',fun.ymax = max, fun.ymin = min, .alpha = 0.05, alpha = 0.5, col="orange",fill="orange")+
		  stat_summary(fun.y=median, geom="line",color="red")+
		  ylim(0,max(s$norm))+
		  theme(text = element_text(size=20))+
		  xlab("Branch min subtree size")+
		  ylab("Ecart")+
		  ggtitle("!{size} Taxa")
	dev.off()
	'''
}

plots.subscribe{
	file -> file.copyTo(resdir.resolve(file.name))
}
