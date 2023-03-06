#purpose is to take a gene ID (RefSeq or such) and create a matrix of gene's N +/- neigbours.
#this can be fed into R or matrix analysis for visualisation
#run in Renv

from os import listdir, getcwd
from os.path import isfile, join, splitext, exists
import pandas as pd
import subprocess
import sys
import pathlib
import argparse



parser = argparse.ArgumentParser(description='Take nucleotide/amino acid sequence, search files for homologues, build a phylo tree')
parser.add_argument('-t','--target_gene', help="A string to look for in attributes(aka 9th aka last) column of GFF file. Ideally gene ID", required=True)
parser.add_argument('-n', '--number_of_genes', default=3, help="Number of genes (on each side of target) to use, default=3", required=False)
parser.add_argument('-d', '--gff_dir', help="The directory with GFF files to search", required=True)
parser.add_argument('-f', '--min_freq', default=3, help="Minimum number of times the gene combination needs to occur to be shown. Prevents excessively large figure. Default=3", required=False)
parser.add_argument('-o', '--output', default="GeneNeighbourhoodDiagram.jpg", help="Output file name. Default=GeneNeighbourhoodDiagram.jpg in current directory", required=False)

args = parser.parse_args()

wd=getcwd()
GffsToExclude=[]
targetGeneID=args.target_gene
neighbourhoodSize=args.number_of_genes
gffDir=args.gff_dir
minComboFrequency=args.min_freq
diagramFileName=args.output

if not exists(gffDir):
    print(f'Could not find directory {gffDir}')
    sys.exit()
files = [f for f in listdir(gffDir) if isfile(join(gffDir, f)) and splitext(f)[1]==".gff"]
if len(files)==0:
    print(f'Could not find files ending in .gff in directory {gffDir}')
    sys.exit()

selectedLines={}

def getGeneLabel(gffLine):
    values=gffLine.split("\t")
    geneDesc=values[8].split(";")
    label="Unknown"
    for i in range(0,len(geneDesc)):
        if "RefSeq:" in geneDesc[i]:
            gene=geneDesc[i].split("RefSeq:")[1]
            label=gene
            #break
        elif "name=" in geneDesc[i]:
            label=geneDesc[i].replace("name=","")
            break
        elif "product=" in geneDesc[i]:
            label=geneDesc[i].replace("product=","")
            break
        elif "note=" in geneDesc[i]:
            label=geneDesc[i].replace("note=","")
            #break
    return label

#not strictly necessary because the logic will assign them, but lint odesn't like unbound variables
targetGenChr=-1
targetGeneDir="-"
#not strictly necessary because the logic will assign them, but lint odesn't like unbound variables

debugCounter=0
for gff in files:
    if gff.replace(".gff","") not in GffsToExclude and gff not in GffsToExclude:
        breakCounter=0 #this is required because the number of lines after the target gene is unknown
        foundGene=False
        lines=[]
        AllLines=[]
        with open(gffDir+gff) as file:
            for line in file:
                if line[0]!="#": #skip the header
                    if line.split("\t")[2].strip()!="CDS":
                        continue
                    AllLines.append(line.strip())
                    if targetGeneID in line:
                        #duplicate diagram markers
                        foundGene=True
                        targetGeneData=line.strip().split("\t")                    
                        targetGenChr=targetGeneData[0]
                        targetGeneDir=targetGeneData[6] #to normalise so that all - target genes are reversed in direction                    
                    if foundGene:
                        breakCounter=breakCounter+1
                        if breakCounter==neighbourhoodSize:
                            lines=lines+AllLines[max(0,len(AllLines)-neighbourhoodSize*2):len(AllLines)]
                            #break #the rest of the file still needs to be parsed because it may have
                            #multiple cases of specific gene
                            foundGene=False
                            breakCounter=0

        #add file data to dictionary
        if foundGene or len(lines)>0:
            selectedLines[gff.replace(".gff","")]=[]
            for line in lines:
                #check that chromosome is the same as the target gene
                neighbourData=line.split("\t")
                if neighbourData[0]==targetGenChr:
                    if targetGeneDir=="-":
                        selectedLines[gff.replace(".gff","")].insert(0, line)
                    else:
                        selectedLines[gff.replace(".gff","")].append(line)
            debugCounter=debugCounter+1
 
df = pd.DataFrame(index=[list(selectedLines.keys())], columns=[("gene"+str(f)) for f in range(-neighbourhoodSize,neighbourhoodSize+1)])
if df.shape[0]==0:
    print("Not found any gff lines with target sequence")
    sys.exit()
for name in selectedLines.keys():
    #name=file name, value=gff lines
    #find position of the target gene in the list of gff lines
    targetGenePos=-1
    for i in range(0,len(selectedLines[name])):
        if targetGeneID in selectedLines[name][i]:
            targetGenePos=i
            break

    lineIndex=0
    for lineIndex in range(0,len(selectedLines[name])):
        line=selectedLines[name][lineIndex]
        #split annotation and get RefSeq, product, or place unknown in that order
        df.at[name,"gene"+str(lineIndex-targetGenePos)]=getGeneLabel(line)

df.to_csv(wd+"/output.tsv", sep="\t") 

#create gggenes input matrix
#the matrix shape is rows are individual genes
#columns are start, end, sample, geneID
#start of KPC is the coordinate 0
gggenesInput=open(wd+"/gggenesInput.tsv","w")
gggenesInput.write("Sample\tgeneID\tStart\tEnd\n")

#not strictly necessary because the logic will assign them, but lint odesn't like unbound variables
flip=True
pointZero=-1
#not strictly necessary because the logic will assign them, but lint odesn't like unbound variables
targetGenChr=-1

for name in selectedLines.keys():
    #name=file name, value=gff lines
    #find position of the target gene in the list of gff lines
    maxPosition=0
    for i in range(0,len(selectedLines[name])):
        values=selectedLines[name][i].split("\t")
        maxPosition=max(maxPosition, int(values[3]),int(values[4]))
        if targetGeneID in selectedLines[name][i]:
            direction=values[6]
            if direction=="+":
                pointZero=int(values[3])
                flip=False
            else:
                pointZero=int(values[4])
                flip=True
            break
    if flip:
        pointZero=maxPosition-pointZero
    for line in selectedLines[name]:
        values=line.split("\t")
        if flip:
            #flip the whole sequence and reverse direction of gene
            values[3]=maxPosition-int(values[3])
            values[4]=maxPosition-int(values[4])

            if values[6]=="+":
                #pass
                gggenesInput.write(name+"\t"+getGeneLabel(line)+"\t"+str(pointZero-int(values[3]))+"\t"+str(pointZero-int(values[4]))+"\n")
            else:
                #pass
                gggenesInput.write(name+"\t"+getGeneLabel(line)+"\t"+str(pointZero-int(values[4]))+"\t"+str(pointZero-int(values[3]))+"\n")
        else:
            pass
            #only adjust start and end coordiantes
            if values[6]=="+":
                gggenesInput.write(name+"\t"+getGeneLabel(line)+"\t"+str(pointZero-int(values[3]))+"\t"+str(pointZero-int(values[4]))+"\n")
            else:
                gggenesInput.write(name+"\t"+getGeneLabel(line)+"\t"+str(pointZero-int(values[4]))+"\t"+str(pointZero-int(values[3]))+"\n")
gggenesInput.close()

subprocess.call (f'Rscript {pathlib.Path(__file__).parent.resolve()}/genesDiagram.R {wd} {diagramFileName} {str(minComboFrequency)}', shell=True)

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread(wd+"/"+diagramFileName)
imgplot = plt.imshow(img)
plt.show()