library(ggplot2)
library(gggenes)
library(dbplyr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("No working directory supplied, check R call in python", call.=FALSE)
} else if (length(args)==3) {
  wd=args[1]
  diagramFile=args[2]
  minComboFrequency=args[3]
}

genes=read.delim(paste0(wd,"/gggenesInput.tsv"),sep="\t")
rejected_assemblies=c("3042","5145","Kp_4240","Kp4151","3766","5170","Kp4189","Kp4173","3763","Kp4179","Kp931")
index=which(genes$Sample %in% rejected_assemblies)
if (length(index)>0){genes=genes[-index,]}

genes$merged=paste(genes$geneID,genes$Start,genes$End)
uniq=unique(genes[,2:ncol(genes)])

#identify unique groups of genes+positions to simplify chart
samples <- list() 
for (sample in unique(genes$Sample)){
  index=which(genes$Sample==sample)
  #the order of uniq$merged %in% genes$merged[index] is very important, reverse order gives non-sense results
  samples[[sample]]=which(uniq$merged %in% genes$merged[index], arr.ind = TRUE)
}


uniqueCombos=unique(samples)
comboCounts=integer(length(uniqueCombos))
for (i in seq(1,length(uniqueCombos))){
  for (sample in samples){
    if (all.equal(sample, uniqueCombos[[i]])==TRUE){
      comboCounts[i]=comboCounts[i]+1
    }
  }
}


figureData=data.frame(Sample=character(0), geneID=character(0), Start=numeric(0), End=numeric(0))
for (i in seq(1,length(uniqueCombos))){
  if (comboCounts[i]>=minComboFrequency){
    newrows=uniq[uniqueCombos[[i]],c("geneID", "Start","End", "merged")]
    #newrows$Sample=paste0("Combo",i,"(",comboCounts[i],")")
    newrows$Sample=paste0("Profile ",i,"(",comboCounts[i],")")
    figureData=rbind(figureData,newrows)
  }
}

figureData$ExtraLabel=""
figureData$geneID=as.character(figureData$geneID) #for some reason, in some R versions this column is made into factor which causes problems
uniqueGenes=unique(figureData[,c("geneID","merged")])
for (i in seq(1,length(uniqueGenes$geneID))){
  index=which(figureData$merged==uniqueGenes$merged[i],arr.ind=TRUE)
  figureData$geneID[index]=paste0(uniqueGenes$geneID[i],"(",as.roman(i),")")
  figureData$ExtraLabel[index]=i
}


p=ggplot(figureData, aes(xmin = Start, xmax = End, y = Sample, fill=geneID, label = as.roman(ExtraLabel))) +
  geom_gene_arrow() + theme(legend.position="bottom", legend.text=element_text(size=8), panel.background = element_rect(fill = 'white', colour = 'white'))+
  geom_gene_label()+
  labs(y="Genes profile", x = "bps")


ggsave(paste0(wd,'/',diagramFile),plot=p, width=40, units=c("cm"), limitsize=FALSE)
#ggsave(paste0(wd,'/','GeneNeighbourhoodDiagram.pdf'),plot=p, width=40, units=c("cm"), limitsize=FALSE)


print("Done with R")