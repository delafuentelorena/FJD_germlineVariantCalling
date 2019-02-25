### CNV analysis 
### Author: Lorena de la Fuente
### Date: 14/01/19


rm(list=ls())

#start clock
ptm <- proc.time()
print(ptm)


#**************************#
### package requirements ###
#**************************#

list.of.packages <- c("ExomeDepth","optparse","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ExomeDepth)
library(GenomicRanges)
library(optparse)
library(stringr)

data(exons.hg19)


# #######################
# opt=NULL
# opt$dir="/home/lorena/floridaCluster/outputAnalysis/test1/applied_bqsr_data"
# opt$dir="/home/lorena/floridaCluster/data/CNV_test_allCancers/applied_bqsr_data"

# opt$outputdir="//home/lorena/floridaCluster/outputAnalysis/test1"
# opt$outputdir="/home/lorena/floridaCluster/data/CNV_test_allCancers"

# opt$bed="/home/lorena/floridaCluster/data/CANCER/18-4175_CA-33/TruSight_Cancer_TargetedRegions_v1.0Used.bed"

# opt$samples="/home/lorena/floridaCluster/data/CNV_test_allCancers/samples.txt"

# #####################


#**************#
#   Arguments  #
#**************#

option_list=list(
  make_option(c("-d","--dir"),type="character", help="Directory with input bam files to analyse."),
  make_option(c('-o','--outputdir'),type="character", help="Output directory."),
  make_option(c('-n','--name'),type="character", help="Output directory."),
  make_option(c('-b','--bed'),type="character",default="genome", help="Bed file with genes found in panel"),
  make_option(c('-s','--samples'), type="character", help="File with negative samples"))


opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args
cnv_output=paste(opt$outputdir,"/copy_number_variation_data",sep="")

sink(paste(opt$outputdir,"/software_",opt$name,".txt",sep=""),append=TRUE)
print("R SESSION INFO (ExomeDepth):")
sessionInfo()
sink() 


#**********************#
# Samples to call CNVs #
#**********************#

if(opt$samples!="all"){
  samples = read.delim2(opt$samples, sep=",", stringsAsFactors = F, header=F)[,2]
}


#*************#
#   Regions   #
#*************#

if(opt$bed=="genome"){
        myexons = exons.hg19
        myexons[,1] = paste("chr",myexons[,1] ,sep = "") 
        }else{
        myexons = read.delim(file=opt$bed, header = F, stringsAsFactors = F)
        colnames(myexons) = colnames(exons.hg19)
        myexons = myexons[-which(myexons$chromosome %in% c("chrX","chrY", "X","Y")),]
        }


cat('\nTotal number of testing regions: ', nrow(myexons), '\n', sep = "")

#reference.fasta <- '/home/lorena/floridaCluster/broad_bundle_hg19_v/ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta'


#******************************#
#   List of sample bam files   #
#******************************#

setwd(opt$dir)
mybams <- Sys.glob('*.bam')
mybais <- Sys.glob('*.bai')



#***********************************#
#   counts dataframe for all BAMs   #
#***********************************#

cat('\nParsing bam files to get counts...\n')

my.counts <- getBamCounts(bed.frame = myexons, 
                          bam.files = mybams,
                          index.files = mybais,
                          min.mapq = 20,
                          include.chr = FALSE)
                          #referenceFasta = reference.fasta) # fasta if we want to compute GC %

setwd(cnv_output)

save.image(file = paste(opt$name,"_counts_image.RData", sep=''))

ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')

ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space),pattern = 'chr',replacement = '')





#******************#
#   count matrix   #
#******************#


ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*.bam')])
rownames(ExomeCount.mat) = ExomeCount.dafr$names
nsamples <- ncol(ExomeCount.mat)
#colnames(ExomeCount.mat) = sapply(colnames(my.counts), function(x) sub(".*applied_bqsr_data_ *(.*?) *.bam*", "\\1", x))
colnames(ExomeCount.mat) = sapply(colnames(my.counts), function(x) sub(".bam", "",strsplit(x, split = "_",  fixed = T)[[1]][1]))

if(opt$samples!="all" && !all(samples %in% colnames(ExomeCount.mat))){
  cat('\nSome of the selected samples can not be found...\n')
  }
if(opt$samples=="all"){
  samples=colnames(ExomeCount.mat) 
}

cat('\nSamples to analyse: ', paste(samples, collapse = ","), '\n', sep = "")
cat('\nAll detected samples: ',paste(colnames(ExomeCount.mat), collapse=','),'\n', sep = "")



#******************#
#   CNV calling    #
#******************#

CNVcalls_sorted_total.toAnnotate = NULL
CNVcalls_sorted_total = NULL

for (i in 1:nsamples) {
  
  samplename = colnames(ExomeCount.mat)[i]
  
  if(samplename %in% samples){
  
    cat("\nAnalysing Sample ",samplename,"...\n",sep = "")
    cat("..................................\n")
    

    #### Create the aggregate reference set for this sample
    
    cat('Read number distribution...\n')
    print(summary(ExomeCount.mat[,i]))
  
    my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                       reference.counts = ExomeCount.mat[,-i,drop=F],
                                       bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                       n.bins.reduced = 10000)
  
    # PCA per sample
    # myfactors=data.frame(row.names = colnames(ExomeCount.mat), sample=colnames(ExomeCount.mat)=="17-2983")
    # myfactors[which(rownames(myfactors)%in%my.choice$reference.choice),"sample"]<-"reference"
    # myfactors[which(rownames(myfactors)==samplename),"sample"]<-"test"
    # mydata2 = readData(ExomeCount.mat,factors = myfactors)
    # myPCA = dat(mydata2, type = "PCA")
    # explo.plot(myPCA, factor = "sample")
    
    cat('Reference set choice: ', paste(my.choice$reference.choice, collapse = ","),"\n",sep = "")
    
    my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                                      MAR = 1,
                                      FUN = sum)
    
    cat('Now creating the ExomeDepth object...\n')
    
    all.exons <- new('ExomeDepth',
                     test = ExomeCount.mat[,i],
                     reference = my.reference.selected,
                     formula = 'cbind(test, reference) ~ 1')

    cat('Calling CNVs...\n')
    
    
    ################ Now call the CNVs
    
    all.exons <- CallCNVs(x = all.exons,
                          transition.probability = 10^-4,
                          chromosome = ExomeCount.dafr$chromosome,
                          start = ExomeCount.dafr$start,
                          end = ExomeCount.dafr$end,
                          name = ExomeCount.dafr$names)
    
    cat("Detected CNVs:",nrow(all.exons@CNV.calls), "\n", sep=" ")
    
    if (nrow(all.exons@CNV.calls) > 0) {
    
      # Adding name found in panel
      
      rangesMyAnnotation = GenomicRanges::GRanges(seqnames = all.exons@annotations$chromosome,
                                                  IRanges::IRanges(start=all.exons@annotations$start,end=all.exons@annotations$end),
                                                  names = all.exons@annotations$name)
      all.exons <- AnnotateExtra(x = all.exons,
                                 reference.annotation = rangesMyAnnotation,
                                 min.overlap = 0.0000001,
                                 column.name =
                                   'annotation')
      
      # Writting file for VEP annotation
      
      CNVcalls_sorted = all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),]
      # Output in tab format ready for VEP
      #forvep.file=paste('CNV_', samplename, '_ready_toAnnotate', '.txt', sep = '')
      
      CNVcalls_sorted[CNVcalls_sorted$type=="duplication","CNV_type"]<-"DUP"
      CNVcalls_sorted[CNVcalls_sorted$type=="deletion","CNV_type"]<-"DEL"
      
      #write.table(x = CNVcalls_sorted[,c("chromosome","start","end","CNV_type"),], file=forvep.file, quote=F, col.names = F, row.names=F, sep="\t")
      
      #output.file=paste('CNV_', samplename, '.txt', sep = '')
      #write.table(file = output.file, x = CNVcalls_sorted[,-13], row.names = FALSE, sep="\t", quote=F)
    
      #CNVcalls_sorted_total = rbind( CNVcalls_sorted_total, data.frame(samplename,CNVcalls_sorted, cor(my.test, my.ref.counts), length( my.choice[[1]])))
      
      CNVcalls_sorted_total = rbind( CNVcalls_sorted_total, data.frame(samplename,CNVcalls_sorted))
      CNVcalls_sorted_total.toAnnotate = unique(rbind(CNVcalls_sorted_total.toAnnotate, CNVcalls_sorted[,c("chromosome","start","end","CNV_type"),]))
  
  
    }
  }

}

#save.image(file = paste(cnv_output,"/",ols pt$name,"_image.RData", sep=''))

output.file = paste('CNV_results_', opt$name, '_exomedepth.txt', sep = '')
forvep.file = paste('CNV_results_', opt$name,  "_toAnnotate", '.txt', sep = '')

cat("\nOUTPUT FILES:\n")
print(output.file)
print(forvep.file)

write.table(x = CNVcalls_sorted_total, file = output.file, row.names = FALSE, sep="\t", quote=F)
write.table(x = CNVcalls_sorted_total.toAnnotate, file=forvep.file, quote=F, col.names = F, row.names=F, sep="\t")

cat("\nCNV DETECTION FINISHED\n")

#stop clock
finish <- proc.time() - ptm
print(finish)
summary(finish)
