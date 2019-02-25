#!/bin/sh

##################################################
##  FJD pipeline - Copy Number Variants Calling ##
##################################################

# REQUIRED MODULES:

VEP_CACHE='/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/.vep'
VEP_FASTA="${VEP_CACHE}/homo_sapiens/93_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
VEP="/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/ensembl-vep/vep"


# CNV input/output directories
MDAP=$1
HCGVCFD=$2
CNV="${MDAP}/copy_number_variation_data"
SAMPLEFILE=$3
run=$4
threads=$5
panel=$6
window=$7
utilitiesPath=$8


softwareFile="${MDAP}/software_${run}.txt"
echo "COPY NUMBER VARIANT CALLING:" >> ${softwareFile}
module list 2>> ${softwareFile}



#####################
## COMPUTING CNVs ###
#####################


echo "mkdir ${MDAP}/copy_number_variation_data"
mkdir $CNV


# modify window size

if [ "$window" != "no" ] &&  [ "$panel" != "genome" ]; then

	panel_out=$(basename "$panel" .bed)_${window}bpwindow.bed
	python $utilitiesPath/CNV_windowSize.py $panel $window $panel_out
	panel=$panel_out
	echo $panel 

fi

# calling CNVs

Rscript $utilitiesPath/exomeDepth.R -d $HCGVCFD -o $MDAP -b $panel -n $run -s $SAMPLEFILE







#####################
## CNVs ANNOTATION ##
#####################


perl $VEP -i $CNV/CNV_results_${run}_toAnnotate.txt \
-o $CNV/CNV_results_VEP_${run}.txt \
--everything --cache --offline --dir $VEP_CACHE --v --assembly GRCh37 \
--fasta $VEP_FASTA --force_overwrite  --tab --refseq



echo -e "\nVEP annotation in $CNV/CNV_results_VEP_${run}.txt"



##############################
## FINAL CNV ANNOTATED FILE ##
##############################


echo -e "\nMerging exomeDepth and VEP output files..."


Rscript $utilitiesPath/vep_processing.R -c $CNV/CNV_results_${run}_exomedepth.txt -v $CNV/CNV_results_VEP_${run}.txt -o $CNV/FINAL_CNV_results_${run}.txt


echo -e "\nFINAL CNV RESULTS in $CNV/FINAL_CNV_results_${run}.txt"
