### VEP file processing 
### Author: Lorena de la Fuente
### Date: 06/02/19


rm(list=ls())

#start clock
ptm <- proc.time()
print(ptm)

library(optparse)
library(stringr)



#**************#
#   Arguments  #
#**************#

option_list=list(
  make_option(c("-c","--cnv"),type="character", help="CNV file per sample"),
  make_option(c('-v','--vep'),type="character", help="Vep text file"),
  make_option(c('-o','--output'),type="character", help="Output file."))

opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args

print(opt$output)

##########

# opt=NULL
# opt$vep="/home/lorena/Documents/TransBioNet/CANCER_TEST/resultsFromBam/copy_number_variation_data/CNV_results_VEP_cftbam_2019_02_06_12_52_20.txt"
# opt$cnv="/home/lorena/Documents/TransBioNet/CANCER_TEST/resultsFromBam/copy_number_variation_data/CNV_results_cftbam_2019_02_06_12_52_20_exomedepth.txt"


##########


#******************#
#   Reading files  #
#******************#

tmp_vcf<-readLines(opt$vep)
header = tmp_vcf[grep(pattern = "^#Upload", tmp_vcf)]
vcf_names<-unlist(strsplit(header,"\t"))
vep = read.table(opt$vep, header=F, quote="\"", stringsAsFactors = F, comment.char = )
colnames(vep) = vcf_names
vepColumns = c("#Uploaded_variation", "Location","Allele","Gene","SYMBOL","BIOTYPE",
               "Feature" ,   "CANONICAL",       "Feature_type" ,   
               "Consequence", "IMPACT",  "Existing_variation" ,"CANONICAL" , "SWISSPROT" ,"cDNA_position") 

vep = vep[,vepColumns]
cnv = read.delim2(opt$cnv, sep="\t", stringsAsFactors = F, header=T)


#******************#
#   Merging files  #
#******************#

cnv$label = paste(cnv$chromosome, ":", cnv$start, "-", cnv$end, ":", cnv$type, sep="")
vep$label = paste(vep$Location, ":", vep$Allele, sep="")

#venn(list(cnv$label, vep$label), cexil =   2)

merged = merge(cnv,vep, by="label", all=T)

merged2 = merged[order(merged$samplename), -which(colnames(merged) %in% c("Allele","id","CNV_type","Location", "#Uploaded_variation")) ]



#******************#
#   Writing file   #
#******************#

print(opt$output)
write.table(x = merged2, file = opt$output, row.names = FALSE, sep="\t", quote=F)



