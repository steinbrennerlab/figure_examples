#source("https://bioconductor.org/biocLite.R")
#biocLite("EBImage")
#biocLite("treeio")
#biocLite("ggtree")
#install.packages("phytools")
#install.packages("optparse")
#install.packages("tidyselect")
#install.packages("tidyselect")
#install.packages("labeling")
#install.packages("devtools")
#devtools::install_github("GuangchuangYu/treeio")

#ggtree script
#see ggtree.Rmd for documentation and a test
#To test locally, download example data and call "Rscript ggtree.R -e tree_input.nwk -o output"

library("ggplot2")
library("treeio")
library("phytools") # for sims and ASRs
library("EBImage") # for images
library("ggtree")
library("optparse")

option_list <- list( 
  make_option(c("-e", "--entry"), action="store", type="character",
              help="entry"),
  make_option(c("-o", "--output"), action="store", type="character", 
              help="name for file writing"),
  make_option(c("-n", "--node"), action="store", default=0, type="integer", 
              help="number of total sequences!  used to compute font size")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir<-paste("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/")
tree <- read.tree(paste(dir,opt$entry,sep=""))
#tree <- read.tree(paste(dir,"tree_input.nwk",sep=""))

dd <- read.table("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/attribute_species.txt", sep="\t", header = TRUE, stringsAsFactor=F)

counts_file <- read.table("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/log2FoldChanges.txt", sep="\t", row.names = 1, header = TRUE, stringsAsFactor=F)

#Set upper and lower limits, used to set colors!
upper <- 5
lower <- -5

#node<-opt$node
node<-0
###
if (node>0) {
  nodenum <- opt$node
  message("Node input detected.  Subtree based on node:")
  message(opt$node)
  tree <- tree_subset(tree, nodenum, levels_back = 1) #Right now a bug with levels_back=0 is preventing me from specifying the node ITSELF
}

message(tree)

q <- ggtree(tree, size=0.1) #size specifies line size
###OPTIONAL: adds hmm coding to ggtree object q
###OPTIONAL: takes hmm_coding and adds to a separate dataframe dd2
dd2 <- read.table("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/attribute_hmm.txt", sep="\t", header = TRUE, stringsAsFactor=F)
q <- q %<+% dd2

node_count <- length(tree$tip.label)
print("The number of nodes is:")
print(node_count)

size <-  3.63 - (0.484*log(node_count)) #computes appropriate font size for tree based on good sizes for 10, 30, and 100
#size <- 2.5 #use a default font size instead
size2 <- (size/2)
print("The tip label font size is")
message(size)
print("The node label font size is")
message(size2)

q <- q %<+% dd + geom_tiplab(size=size,offset=0.05,aes(color=species)) + geom_tippoint(aes(size=size2,shape=hmm)) + scale_size_identity() #you need scale_size_identity! https://groups.google.com/forum/#!topic/bioc-ggtree/XHaq9Sk3b00

figure <- gheatmap(q,counts_file, offset = 1.5, width=0.6, font.size=1.5, colnames_angle=-45, hjust=0) + 
  geom_tiplab(size=size,offset=0.05,aes(color=species)) +
  scale_fill_gradient2(low="#000099",high="#FF0000",mid="white",limits=c(lower,upper)) +
  #node labels
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=size2) + 
  scale_colour_manual(values=c("black","red","blue","orange","purple","darkgreen","cadetblue","deeppink","darkgoldenrod","brown4","olivedrab2"))
figure


file <- paste(dir,opt$output,".pdf",sep="")
#file <- "C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/ggtree_output.pdf"
message(file)
pdf(file)
figure
dev.off()

#Takes the tree object and converts it to a tidy dataframe using fortify, then reorders it according to the graphical position!
#Apparently fortify might deprecate and switch to the "broom" package for tidying data
tips <- fortify(tree)
tips <- data.frame(tips$label,tips$y,tips$isTip)
tips <- tips[order(tips[,3],tips[,2],decreasing=T),]
#Writes the tips to a csv file.  Name is based on the option -b specified when the script is called
file_csv <- paste(dir,opt$output,".csv",sep="")
#file_csv <- "C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/ggtree_output.csv"
message(file_csv)

for(i in 1:node_count) {
  write(as.matrix(tips)[,1][i],file=file_csv,sep=",",append=T)
}