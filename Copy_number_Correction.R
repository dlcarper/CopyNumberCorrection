library(rstudioapi)
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path ))
print( getwd() )
library(phyloseq)
library(ape)
library(vegan)
library(ampvis)
library(viridis)
library(ggpubr)
library(ggplot2)
library(readxl)

#Import into phyloseq
otu<-read.table("feature-table.txt",sep="\t",header=TRUE)
Taxa<-read.table("Taxa_table1.txt",sep="\t",header=TRUE)

otu2<-merge(otu,Taxa, by="OTU.ID")
copynumber<-read_excel("16S_copyNumber.xlsx")
otu2_withcopy<-merge(otu2,copynumber,by.x="Species",by.y = "Sample.Name",all.x = TRUE)
rownames(otu2_withcopy)<-paste0(otu2_withcopy$OTU.ID)
otu2_withcopy$OTU.ID<-NULL
otu2_withcopy$IMG.Genome.ID<-as.factor(otu2_withcopy$IMG.Genome.ID)
otu2_numeric<-select_if(otu2_withcopy,is.numeric)
otu3<-otu2_numeric/otu2_numeric$`16S.rRNA.Count`
otu3$`16S.rRNA.Count`<-NULL
otumatrix<-as.matrix(otu3)

rownames(Taxa)<-paste0(Taxa$OTU.ID)
Taxa$OTU.ID<-NULL
taxamatrix<-as.matrix(Taxa)
Mapping<-read.table("Test_Sequencing_mapping_file_R.txt",sep="\t",header=TRUE,row.names = 1)
OTU = otu_table(otumatrix, taxa_are_rows = TRUE)
TAX = tax_table(taxamatrix)
tree = read.tree("tree.nwk")
META = sample_data(Mapping)
physeq1 = phyloseq(OTU,TAX,META,tree)
