#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)>4) {
  stop("Usage: Rscript plot_identity.r input_file addline_onY outputPDF [highlight.fa]", call.=FALSE)
} 

dat = read.delim(args[1], header=F, sep="")
hline_y = as.numeric(args[2])

df_hlRegions <- try(read.delim(args[4], header=F, sep=""), silent=TRUE)

#list_colors <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')
list_colors <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')

library(ggplot2)
	pdf(args[3], width=12, height=8)
if (class(df_hlRegions) == "try-error"){
	ggplot(data=dat, aes(x=(V1+V2)/2, y=V3, colour=V4)) + geom_line() + geom_hline(yintercept=hline_y, linetype="dashed") + facet_wrap(~V4,ncol=1) + xlab("Position on alignment") + ylab("Fraction of identity") + theme_bw() + theme(legend.position="bottom", strip.text.x = element_text(size = 13), axis.text.x= element_text(size=12), axis.text.y=element_text(size=13),axis.title.x = element_text(size=13), axis.title.y = element_text(size=13), legend.text=element_text(size=13), legend.title=element_text(size=13), plot.title=element_text(size=13))
} else{
	names(df_hlRegions) <- c("hapID","start","end","dupID")
	df_hlRegions$label_x = (df_hlRegions$start + df_hlRegions$end) / 2
	df_hlRegions$label_y = (max(dat$V3) + min(dat$V3))/2
	df_hlRegions$id = paste(paste(df_hlRegions$hapID, df_hlRegions$start, df_hlRegions$end, sep="_"), df_hlRegions$dupID, sep="\n")
	df_hlRegions$dupID <- factor(df_hlRegions$dupID)
	mapping_colors <- data.frame(type=levels(df_hlRegions$dupID), color=list_colors[1:length(levels(df_hlRegions$dupID))])
	df_hlRegions$color <- mapping_colors[ match(df_hlRegions$dupID, mapping_colors$type), "color"]

	ggplot(data=dat, aes(x=(V1+V2)/2, y=V3, colour=V4)) + geom_line() + 
	geom_rect(data=df_hlRegions, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(dat$V3), ymax=1), fill=rep(df_hlRegions$color,3), alpha=0.3) +
	geom_text(data=df_hlRegions, inherit.aes=FALSE, aes(x=rep(df_hlRegions$label_x,3), y=rep(df_hlRegions$label_y,3), label=rep(df_hlRegions$id,3)), size=2.5, angle=90) + 
	geom_hline(yintercept=hline_y, linetype="dashed") + facet_wrap(~V4,ncol=1) + xlab("Position on alignment") + ylab("Fraction of identity") + theme_bw() + theme(legend.position="bottom", legend.title = element_blank(), strip.text.x = element_text(size = 13), axis.text.x= element_text(size=12), axis.text.y=element_text(size=13),axis.title.x = element_text(size=13), axis.title.y = element_text(size=13), legend.text=element_text(size=13), plot.title=element_text(size=13))  

}
dev.off()
