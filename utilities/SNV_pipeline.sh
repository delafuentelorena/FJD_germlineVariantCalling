#!/bin/sh

################################
### FJD pipeline - Haplotype ###
################################

### FJD PIPELINE ARGUMENTS:

INPUT=$1
MDAP=$2
sample=$3
threads=$4
run=$5
panel=$6
analysis=$7
cvcf=$8
skipmapping=${9}    
HG19=${10}
local=${11}


# REQUIRED MODULES AND DATABASES:

if [ "$local" != "True" ]; then

	module load bwa/0.7.17
	module load samtools/1.3.1
	module load picard/2.18.20
	module load gatk/4.0.8.0
	module load vep/92
	module load bedtools/2.24.0
	module load R
	alias picard='java -jar $HPC_PICARD_DIR/picard.jar'
	alias gatk='java -jar $HPC_GATK_DIR/gatk-package-4.0.8.0-local.jar'


else

	export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
	alias bwa='$SFT/bwa/bwa'
	alias samtools='$SFT/samtools/samtools'
	alias picard='java -jar $SFT/picard/build/libs/picard.jar'
	alias gatk='java -jar $SFT/gatk/build/libs/gatk-package-4.0.6.0-22-g9d9484f-SNAPSHOT-local.jar'
	VEP="${SFT}/variant_effect_predictor/ensembl-vep/vep"
	FILTER_VEP="${SFT}/variant_effect_predictor/ensembl-vep/filter_vep"
	VEP_CACHE='/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/.vep'
	VEP_FASTA="${VEP_CACHE}/homo_sapiens/93_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
	PLUGIN_DIR="${VEP_CACHE}/Plugins"
	PLUGIN_DBS="${VEP_CACHE}/dbs"
	DBNSFP="${PLUGIN_DBS}/dbNSFP_hg19.gz"

fi




if [ "$skipmapping" != "True" ]; then

	foward="${INPUT}/${sample}*R1*.fastq.gz"
	reverse="${INPUT}/${sample}*R2*.fastq.gz"

else

	mv ${INPUT}/${sample}*bam ${INPUT}/${sample}_alignment.bam
	mv ${INPUT}/${sample}*bai ${INPUT}/${sample}_alignment.bai

fi


softwareFile="${MDAP}/software_${run}.txt"
echo "MAPPING / SINGLE NUCLEOTIDE VARIANT CALLING:" >> ${softwareFile}
module list 2>> ${softwareFile}



### FJD PIPELINE VARIABLES:
#·······································

#Mapped_data
MD="${MDAP}/mapped_data"

#Sorted_data
SD="${MDAP}/sorted_data"
	
#Dedupped_data
DD="${MDAP}/dedupped_data"
#Recalibrated_base_quality_scores_data: generates a recalibration table  based on specified covariates, read_group,reported_quiality_score,machine_cycle and nucleotide_context, after GATK.
RBQSRD="${MDAP}/recalibrated_bqsr_data"
#Applied_bqsr_data: applying the reacalibration table to the BAM file to continue the analysis based on the READS best selected by GATK.

if [ "$skipmapping" != "True" ]; then
	ABQSRD="${MDAP}/applied_bqsr_data"
	mkdir $ABQSRD
else
	ABQSRD="${INPUT}"
fi

#Plot_recalibration_data: plots the recalibration differences between the first and the second pass of the recalibration.
#First using the APPLIED_RECALIBRATION_BQSR_DATA BAM file, running again BaseRecalibrator and generating the plots with AnalyzeCovariates.
PRD="${MDAP}/plot_recalibration_data"
#Haplotype_caller_gvcf_data:calling for SNPs and indels via local re-assembly of HAPLOTYPES using HAPLOTYPECALLER by GATK.
HCGVCFD="${MDAP}/haplotype_caller_gvcf_data"
#Perform joint genotyping on one or more samples precalled with Haplotype_caller, if one sample,
#straight after haplotypecaller and if more than one use Combine_gvcfs.
#Genotyped_vcf_data: single or single-multisample GVCF as input, output will be a VCF.
GVCFD="${MDAP}/genotyped_vcf_data"
#HARD_FILTERING, filters for SNPs and filters for INDELs.
#Variant_filtration_vcf_data:hard filtering process to select based on the INFO and FORMAT annotations(QD,MQO,FS,MQ,ReadPosrankSum)
VFVCFD="${MDAP}/variant_filtration_vcf_data"
#Vep_vcf_annotated_data: annotations added to the CSQ tag in INFO columnd to the VCF format, --vcf argument (change name in VCF_out)
#if TSV will come out with a separated tab value for each annotation in a column. --tab argument for TSV format
VEPVCFAD="${MDAP}/vep_vcf_annotated_data"

TMP=$MDAP/tmp_${sample}



#Start pipeline processing.

if [ "$skipmapping" != "True" ]; then



	echo "\n\n\n- INDEXING REFERENCE FILES (BWA) "
	echo "-----------------------------------\n"

	# BWA index
	if [ ! -f $HG19/ucsc.hg19.fasta.sa ]; then
		#Start BWA INDEX.
		echo 'Starts BWA INDEX'
		bwa index $HG19/ucsc.hg19.fasta
		echo sampleFile '\nBWA INDEX DONE' 
	else
		echo 'BWA INDEX already existing. Reference indexing for BWA skipped'
	fi




	# Genome index
	echo 
	if [ ! -f $HG19/ucsc.hg19.fasta.fai ]; then
		#Creating .FAI in HG19.
		echo 'Create .FAI file, using samtools faidx'
		samtools faidx $HG19/ucsc.hg19.fasta 
		echo sampleFile '\nucsc.hg19.fasta.FAI DONE'
	else
		echo 'REFERENCE INDEX already existing. Fasta indexing skipped'


	fi



	# Genome dict
	if [ ! -f $HG19/ucsc.hg19.dict ]; then
		#Creating .DICT in HG19.
		echo 'Create .DICT file, using picardtools CreateSequnceDictionary'
		picard CreateSequenceDictionary R=$HG19/ucsc.hg19.fasta O=$HG19/ucsc.hg19.2.dict
		echo sampleFile 'ucsc.hg19.DICT DONE'
	else
		echo 'REFERENCE DICT already existing. Dictionary generation skipped'

	fi







	echo "\n\n\n- MAPPING READS (BWA) "
	echo "-----------------------------\n"

	#mapping fastq files to reference_genome after BWA INDEX
	#reference genome is: UCSC.HG19.FASTA

	mkdir $MD
	echo 'mkdir mapped_data' 

	echo '\nBuilding the header for '${sample}' ongoing...\n'
	#build header
	header=$(zcat $foward | head -n 1)
	echo $header
	id=$(echo $header | head -n 1 | cut -f 1-4 -d':' | sed 's/@//' | sed 's/:/_/g')
	echo $id
	sm=$(echo $header | head -n 1 | grep -Eo '[ATGCN]+$')
	echo $sm
	echo  "\nThis is how the new header looks\n"
	echo '@RG\tID:'$id'\tSM:'${sample}'\tLB:'$id'_'$sm'\tSM:'$id'_'$m'\tPL:ILLUMINA'
	echo  '\nHEADER --> DONE\n'

	bwa mem -v 3 -t $threads -R '@RG\tID:'$id'\tSM:'${sample}'\tLB:'$id'_'$sm'\tPL:ILLUMINA' \
	$HG19/ucsc.hg19.fasta \
	$foward \
	$reverse > $MD/mapped_${sample}.sam
	echo  '\nBWA MEM '${sample}'  DONE'







	echo "\n\n\n- SORTING BAM (PICARD) "
	echo "-------------------------\n"
	mkdir $SD
	mkdir $TMP
	echo 'mkdir sorted_data'
	#Sorting mapped data.
	echo 'Run picard SortSam for '${sample}''
	picard SortSam TMP_DIR=$TMP I=$MD/mapped_${sample}.sam  \
	         O=$SD/sorted${sample}.bam \
	         SO=coordinate
	echo '\nPicard SortSam for '${sample}'  DONE' 









	echo "\n\n\n- MARKING DUPLICATES (PICARD)"
	echo "---------------------------------\n"

	#Selecting the duplicates reads from the mapped and sorted reads.
	mkdir $DD
	echo 'mkdir dedupped_data'
	#Mark duplicates PICARD
	echo 'Start picard MarkDuplicates '${sample}' '
	picard MarkDuplicates TMP_DIR=$TMP \
	 I=$SD/sorted${sample}.bam \
	 O=$DD/dedupped_${sample}.bam \
	 M=$DD/marked_dup_metrics_${sample}.txt \
	 REMOVE_DUPLICATES=true \
	 AS=SortOrder
	echo  '\n PICARD MarkDuplicates '${sample}' DONE' 







	echo "\n\n\n- Dedupped BAM (PICARD)"
	echo "-------------------------\n"

	#Create a .BAI file to compare original vs removed duplicates reads.
	#	#Indexing the BAM files.
	#	#Generating the .BAI files from the DEDUPPED (markedDuplicates from the original SAM/BAM file).
	echo 'Indexing '${sample}' BAM files'
	picard BuildBamIndex TMP_DIR=$TMP \
	 I=$DD/dedupped_${sample}.bam \
	 O=$DD/dedupped_${sample}.bai
	echo  '\n PICARD BuildBamIndex '${sample}' DONE'







	echo "\n\n\n- BQSR RECALIBRATION TABLE (GATK)"
	echo "-------------------------------------\n"


	#Recalibrating  the reads using base quality score reads.
	mkdir $RBQSRD
	echo 'mkdir recalibrated_bqsr_data' 

	#GATK BaseRecalibration first table
	#BaseRecalibration + table
	echo 'Starts GATK '${sample}' Recalibrator'
	gatk BaseRecalibrator --TMP_DIR=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-I $DD/dedupped_${sample}.bam \
	--known-sites $HG19/dbsnp_138.hg19.vcf \
	--known-sites $HG19/1000G_phase1.indels.hg19.sites.vcf \
	--known-sites $HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-O $RBQSRD/before_recalibrated_bqsr_data_${sample}.recal.table
	#--bqsr 1st_racalibration.table

	echo  '\n GATK BaseRecalibrator '${sample}' DONE' 







	echo "\n\n\n- APPLYING BQSR  (GATK) "
	echo "--------------------------\n"


	#Applying the recalibration table to the bam file to continue the analysis.
	#ApplyBQSR
	echo 'Starts picard  '${sample}' ApplyBQSR'
	gatk ApplyBQSR --TMP_DIR=$TMP \
	-R $HG19/ucsc.hg19.fasta \
	-I $DD/dedupped_${sample}.bam \
	--bqsr $RBQSRD/before_recalibrated_bqsr_data_${sample}.recal.table \
	-O $ABQSRD/${sample}_alignment.bam
	echo  '\nGATK ApplyBQSR '${sample}' DONE'


	echo "\n\n\n- ANALYZE COVARIATES FOR PDF PLOTS. COMPARISON OF THE RECALIBRATION DATA (GATK)"
	echo "----------------------------------------------------------------------------------\n"


	mkdir $PRD
	echo 'mkdir plot_recalibration_data' 
	echo 'Starts GATK Second Recalibration'

	#GATK BaseRecalibration second table for next step AnalyzeCovariates.
	#Generates the second pass table.
	#instead of second table, we use the BAM created by ApplyBQSR to regenrate a new TABLE for plot.

	echo $ABQSRD/${sample}_alignment.bam
	# gatk BaseRecalibrator --TMP_DIR=$TMP \
	# -I $ABQSRD/${sample}_alignment.bam \
	# -R $HG19/ucsc.hg19.fasta \
	# --known-sites $HG19/dbsnp_138.hg19.vcf \
	# --known-sites $HG19/1000G_phase1.indels.hg19.sites.vcf \
	# --known-sites $HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	# -O $PRD/after_recalibrated_bqsr_data_${sample}.recal.table
	# echo 'Second recalibration GATK table --> DONE'

		#Full generating the Plots of the recalibration tables using AnalyzeCovariates and saving a csv copy.
		#Analyze  the tables.

	echo   '\n  Generating Plots, pdf and csv files...'
	# gatk AnalyzeCovariates --TMP_DIR=$TMP \
	# -before $RBQSRD/before_recalibrated_bqsr_data_${sample}.recal.table \
	# -after $PRD/after_recalibrated_bqsr_data_${sample}.recal.table \
	# -csv $PRD/BQSR_${sample}.csv \
	# -plots $PRD/AnalyzeCovariates_bqsr_${sample}.pdf
	# echo '\nPlot and CSV file for ${sample} is DONE!'
	echo sampleFile '\n Plots files generated --> DONE'
	#Obtaining an CSV and PDF file of the comparrisson between first and second pass of the recalibration applied to the bam fi$



	rm $MD/mapped_${sample}.sam
	rm $SD/sorted${sample}.bam
	rm $DD/dedupped_${sample}.bam
	rm $DD/marked_dup_metrics_${sample}.txt
	rm $DD/dedupped_${sample}.bai



fi



### If just mapping step, then exit



if [ "$analysis" = "mapping" ]; then

	rm -r $TMP
	exit 0

fi






echo "\n\n\n- HAPLOTYPECALLER (GATK)"
echo "---------------------------\n"


#Ready to call for Variants.
mkdir $HCGVCFD
echo 'mkdir haplotype_caller_gvcf_data' 

#HaplotypeCaller for each sample for later joint genotyping.
echo sampleFile '\nGATK HaplotypeCallerGVCF for '${sample}' STARTS'
gatk HaplotypeCaller --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-I $ABQSRD/${sample}_alignment.bam \
-ERC GVCF \
-bamout $HCGVCFD/HC_bamout_${sample}.bam \
-O $HCGVCFD/HC_data_${sample}.g.vcf \
-G StandardAnnotation \
-G AS_StandardAnnotation \
-G StandardHCAnnotation \
-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet -A StrandArtifact \
--annotate-with-num-discovered-alleles=true 


# if [ "$panel" != "genome" ]; then
# 	bgzip -c $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf > $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf.gz
# 	bcftools index $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf.gz
# 	bcftools filter -R $panel -o $VFVCFD/pre_filtered_ready_to_annotate_TARGETING_${sample}.vcf $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf.gz 
# 	bcftools stats 	$VFVCFD/pre_filtered_ready_to_annotate_TARGETING_${sample}.vcf > $VFVCFD/pre_filtered_ready_to_annotate_TARGETING_${sample}.stats
# 
#else
#
#fi

echo sampleFile '\nGATK HaplotypeCallerGVCF ERC GVCF for '${sample}' DONE' 

echo $HCGVCFD/HC_data_${sample}.g.vcf >> $HCGVCFD/my_list_of_gvcfs_files_to_combine_$run.list




# if joint vcf exit here


if [ "$cvcf" = "True" ]; then

	exit 0

fi








echo "\n\n\n- GENOTYPECALLER (GATK) "
echo "---------------------------\n"


#GenotypeGVCFs into final VCF

mkdir $GVCFD
echo 'mkdir genotyped_data_vcf'

#GenotypeGVCF into final VCF
echo sampleFile '\nUsing GATK GenotypeGVCFs for final VCF'
gatk GenotypeGVCFs --TMP_DIR=$TMP \
    -R $HG19/ucsc.hg19.fasta \
    -V $HCGVCFD/HC_data_${sample}.g.vcf \
    -G StandardAnnotation \
    -O $GVCFD/genotyped_data_${sample}.vcf 
    
echo sampleFile '\nGATK GenotypeGVCFs for '${sample}' DONE'










echo "\n\n\n- HARD FILTERING (classical way GATK)"
echo "----------------------------------------\n"


#HARD FILTERING
#First step extacting the SNP's
#Second step extracting the INDEL's


mkdir $VFVCFD
echo "mkdir variant_filtration_data_vcf"


### SNPs

#1.Extract the SNP's from the call set.
echo "Extract the SNP's from the call set."

gatk SelectVariants --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $GVCFD/genotyped_data_${sample}.vcf \
--select-type-to-include SNP \
-O $VFVCFD/selected_raw_snp_${sample}.vcf
#Creates the selected_raw_snps vcf containing just the SNP's from the original  callset.


#2.Apply the filters to the SNP's callset.

gatk VariantFiltration --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $VFVCFD/selected_raw_snp_${sample}.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_SNP_filter" \
-O $VFVCFD/filtered_SNP_data_${sample}.vcf



### INDELs

#3. Extract the INDELS from the ORIGINAL call set.
gatk SelectVariants --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $GVCFD/genotyped_data_${sample}.vcf \
--select-type-to-include INDEL \
-O $VFVCFD/selected_raw_indels_${sample}.vcf
#Creates the selected_raw_indels vcf containing just the INDEL's from the original  callset.


#4.Apply the filters to the INDEL's callset.
gatk VariantFiltration --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $VFVCFD/selected_raw_indels_${sample}.vcf \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "my_INDEL_filter" \
-O $VFVCFD/filtered_INDEL_data_${sample}.vcf


# Combine Variants after using SNPS and INDELS filtering into a single file and get it ready for annotation.

gatk MergeVcfs --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-I $VFVCFD/filtered_SNP_data_${sample}.vcf \
-I $VFVCFD/filtered_INDEL_data_${sample}.vcf \
-O $VFVCFD/filtered_INDEL_SNP_data_${sample}.vcf




# echo ············································································································
# echo -e "Filtering QUALITY,chrM,chrUn,my_SNP_filter and my_INDEL_filter out of the VCF file"
# echo ············································································································

# awk 'BEGIN {FS=OFS="\t"}{if ($1 ~/^#/){print $0} else if ($1 ~/^chrM.*/ || $1 ~/^chrUn.*/  || $7 ~/my_SNP_filter/ || $7 ~/my_INDEL_filter/ || $6 < 100) {print ""} else {print $0}}'  $VFVCFD/filtered_INDEL_SNP_data_${sample}.vcf > $VFVCFD/test_${sample}; gawk 'NF > 0' $VFVCFD/test_${sample} >  $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf

# bcftools stats $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf > $VFVCFD/pre_filtered_ready_to_annotate_${sample}.stats


# #### Filtering by panel coordinates

# if [ "$panel" != "genome" ]; then
# 	bgzip -c $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf > $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf.gz
# 	bcftools index $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf.gz
# 	bcftools filter -R $panel -o $VFVCFD/pre_filtered_ready_to_annotate_TARGETING_${sample}.vcf $VFVCFD/pre_filtered_ready_to_annotate_${sample}.vcf.gz 
# 	bcftools stats 	$VFVCFD/pre_filtered_ready_to_annotate_TARGETING_${sample}.vcf > $VFVCFD/pre_filtered_ready_to_annotate_TARGETING_${sample}.stats
# fi







echo "\n\n\n- VARIANT ANNOTATION (VEP - ENSEMBL) "
echo "----------------------------------------\n"



mkdir $VEPVCFAD
echo "mkdir vep_vcf_annotated_data"

VCF_IN="${VFVCFD}/filtered_INDEL_SNP_data_${sample}.vcf"
VCF_OUT="${VEPVCFAD}/vep_qual_filters_${sample}.vcf"

VCF_FILTERED="${VEPVCFAD}/vep_qual_filters_annotated_${sample}.vcf"
VCF_FILTERED_2="${VEPVCFAD}/vep_qual_filters_annotated_can_conseq_${sample}.vcf"
VCF_FILTERED_3="${VEPVCFAD}/vep_qual_filters_annotated_can_conseq_freqs_${sample}.vcf"
VCF_FINAL="${VEPVCFAD}/FinalVariantAnnotation_${sample}.txt"





echo '\nFiltering by quality and chromosome...\n'

perl $FILTER_VEP \
-i $VCF_IN -o $VCF_OUT \
--filter "QUAL > 100 and FILTER = PASS" \
--filter "CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chrX,chrY" --force_overwrite \





echo '\nVEP annotation...\n'

perl $VEP \
--cache --offline --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR --v --fork 16 --assembly GRCh37 --fasta $VEP_FASTA --force_overwrite \
--biotype --regulatory --protein --symbol --allele_number --numbers --domains --uniprot --variant_class \
--canonical \
--sift p --polyphen p --af --max_af \
--format vcf --vcf \
--pubmed \
--plugin dbscSNV,$PLUGIN_DBS/dbscSNV1.1_GRCh37.txt.gz \
--plugin LoFtool,$PLUGIN_DIR/LoFtool_scores.txt \
--plugin ExACpLI,$PLUGIN_DIR/ExACpLI_values.txt \
--plugin dbNSFP,$DBNSFP,gnomAD_exomes_AF,gnomAD_exomes_NFE_AF,gnomAD_genomes_AF,1000Gp3_AF,1000Gp3_EUR_AF,ExAC_AF,ExAC_EAS_AF,ExAC_NFE_AF,ExAC_Adj_AF,rs_dbSNP150,phyloP20way_mammalian,phyloP20way_mammalian_rankscore,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,GERP++_RS,GERP++_RS_rankscore,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,MetaLR_pred,MetaSVM_pred,M-CAP_pred,Interpro_domain,GTEx_V6p_gene,GTEx_V6p_tissue \
-i $VCF_OUT -o $VCF_FILTERED




echo '\nFiltering by population frequencies...\n'

perl $FILTER_VEP \
-i $VCF_FILTERED -o $VCF_FILTERED_2 \
--filter "DP > 10" \
--filter "CANONICAL = YES" \
--filter "Consequence is 3_prime_UTR_variant or Consequence is 5_prime_UTR_variant or Consequence is intron_variant or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is splice_region_variant or Consequence is synonymous_variant or Consequence is missense_variant or Consequence is inframe_deletion or Consequence is inframe_insertion or Consequence is stop_gained or Consequence is frameshift_variant or Consequence is coding_sequence_variant" --force_overwrite \

perl $FILTER_VEP \
-i $VCF_FILTERED_2 -o $VCF_FILTERED_3 \
--filter "AF < 0.01 or not AF" \
--filter "(ExAC_EAS_AF < 0.01 or not ExAC_EAS_AF) and (1000Gp3_AF < 0.01 or not 1000Gp3_AF) and (1000Gp3_EUR_AF < 0.01 or not 1000Gp3_EUR_AF ) and (gnomAD_exomes_AF < 0.01 or not gnomAD_exomes_AF) and (ExAC_AF < 0.01 or not ExAC_AF) and (ExAC_Adj_AF < 0.01 or not ExAC_Adj_AF) and (ExAC_NFE_AF < 0.01 or not ExAC_NFE_AF) and (gnomAD_exomes_NFE_AF < 0.01 or not gnomAD_exomes_NFE_AF) and (gnomAD_genomes_AF < 0.01 or not gnomAD_genomes_AF)" --force_overwrite \






echo "\n\n\n- OUTPUT PROCESSING "
echo "-------------------------\n"


utilitiesPath=$(dirname "$0")


echo '\nFrom VEP annotated VCF to txt file...\n'


python $utilitiesPath/vep2tsvlorena.py $VCF_FILTERED_3 \
-o $VCF_FINAL -t -g -f GT -i


# Remove temporary files.

echo '\n'
rm -r $TMP




# if [ homo = "True" ]; then


# echo "\n\n\n- HOMOZYGOSITY (PLINK) "
# echo "-------------------------\n"


# 	plink --vcf -recode -out mycancerplink2
# 	../plink --file hapmap1 --homozyg  --homozyg-snp 50 --homozyg-window-het 0
# 	vcftools --vcf genotyped_data_allsamples.vcf --out mycancerplink2 --plink
# 	../plink --bfile hola-temporary --homozyg  --homozyg-snp 50  --homozyg-window-het 0

# 	# Plink default: 
# 	--homozyg-kb 1,000
# 	--homozyg-snp 100

# 	# Plink improvemements (minimizes the trade-off between the exclusion of non-autozygous ROHs and the cost of missing shorter autozygous ROHs)
# 	# - 1. remove SNPs with MAF < 0.05
# 	# - 2. LD prunning in two levels: ligth to moderate outperforms ROH detection -> --indep 
# 	# Parameters of light prunning: r2>0.9 in a 50 SNP window (VIF > 10)
# 	# Parameters of moderate prunning: 50 SNP window that had r2>0.5 (VIF > 2)

# 	# - type1 and type2 error rates are computed at SNP level.

# 	# You may also want to try 'bcftools roh', which uses a HMM-based detection method. 
# 	#(We'll include a basic port of that command in PLINK 2.0 if there is sufficient interest.)

# 	#(For more accurate detection of smaller segments, one might consider approaches that also take
# 	# population parameters such as allele frequency and recombination rate into account, 
# 	# in a HMM approach for example: but for now, PLINK only supports this basic detection of long, homozygous segments). 

# 	#   The originality of this algorithm consists in its capabil-ity to incorporate distances between consecutive SNPs to discriminate between the homozygosity and the hetero-zygosity states. This feature of  H3M2  offers advantage in the detection of small-sized (<500 kb) and medium-sized (between 500 and 1,500 kb) ROH  [18] . 
# 	# Differently,  H3M2and PLINK show similar performances in the WES-based detection of long ROH (>1,500 kb)

# fi

