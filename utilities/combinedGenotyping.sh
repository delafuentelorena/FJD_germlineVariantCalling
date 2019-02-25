#!/bin/sh

############################################################
### FJD pipeline - Combined Genotyping and Trio analysis ###
############################################################


### FJD PIPELINE ARGUMENTS:

MDAP=$1
run=$2
threads=$3
HG19=$4
ped=$5
local=$6
echo $local
echo $1
echo $2
echo $3
echo $4
echo $5




# REQUIRED MODULES AND DATABASES:

if [ "$local" != "True" ]; then
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
	alias gatk3='java -jar /mnt/genetica3/GeneticaPipeDB_updated/gatk-3.8/GenomeAnalysisTK.jar'

	VEP="${SFT}/variant_effect_predictor/ensembl-vep/vep"
	FILTER_VEP="${SFT}/variant_effect_predictor/ensembl-vep/filter_vep"
	VEP_CACHE='/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/.vep'
	VEP_FASTA="${VEP_CACHE}/homo_sapiens/93_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
	PLUGIN_DIR="${VEP_CACHE}/Plugins"
	PLUGIN_DBS="${VEP_CACHE}/dbs"
	DBNSFP="${PLUGIN_DBS}/dbNSFP_hg19.gz"

fi





### VARIABLES:

#Haplotype_caller_gvcf_data:calling for SNPs and indels via local re-assembly of HAPLOTYPES using HAPLOTYPECALLER by GATK.
HCGVCFD="${MDAP}/haplotype_caller_gvcf_data"

#If analyzing more than one sample or an family trio, use COMBINE_GVCFS to  combine them into a single GVCF.
CGVCFD="${MDAP}/combined_gvcf_data"

#Genotyped_vcf_data: single or single-multisample GVCF as input, output will be a VCF.
GVCFD="${MDAP}/genotyped_vcf_data"
	
#Variant_filtration_vcf_data:hard filtering process to select based on the INFO and FORMAT annotations(QD,MQO,FS,MQ,ReadPosrankSum)
VFVCFD="${MDAP}/variant_filtration_vcf_data"

#Vep_vcf_annotated_data: annotations added to the CSQ tag in INFO columnd to the VCF format.
VEPVCFAD="${MDAP}/vep_vcf_annotated_data"


TMP=$MDAP/tmp_$run
mkdir $TMP








### RUN:


softwareFile="${MDAP}/software_${run}.txt"
echo "Combined genotyping:\n" >> ${softwareFile}
module list 2>> ${softwareFile}

echo "HAPLOTYPE DIR: $MDAP"
echo "Run name: $run"







echo "\n\n\n- COMBINED GVCFS (GATK) "
echo "-----------------------------\n"

echo 'mkdir combined_gvcf_data'
mkdir $CGVCFD
echo '\nUsing GATK COMBINEGVCFs for merging GVCFs'

gatk CombineGVCFs --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
--variant $HCGVCFD/my_list_of_gvcfs_files_to_combine_$run.list \
-O $CGVCFD/combined.g.vcf

echo '\nGATK COMBINEGVCFs DONE'











echo "\n\n\n- JOINT GENOTYPING (GATK) "
echo "-----------------------------\n"


mkdir $GVCFD
echo 'mkdir genotyped_data_vcf'
echo '\nUsing GATK GenotypeGVCFs for final VCF'

gatk GenotypeGVCFs --TMP_DIR=$TMP  \
-R $HG19/ucsc.hg19.fasta \
-V $CGVCFD/combined.g.vcf \
-O $GVCFD/genotyped_data_$run.vcf \
-ped $ped

echo '\nGATK GenotypeGVCFs DONE' 












echo "\n\n\n- HARD FILTERING (GATK) "
echo "-----------------------------\n"

echo "\tHard filtering if less than 30 samples and doing it in the classical way, selecting variants\n"
echo "mkdir variantfiltration_data_vcf"
mkdir $VFVCFD


### SNPs
echo "Extract and filter SNP's from the call set."

	#1.Extract the SNP's from the call set.
gatk SelectVariants --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $GVCFD/genotyped_data_$run.vcf \
--select-type-to-include SNP \
-O $VFVCFD/selected_raw_snp_$run.vcf

	#2.Apply the filters to the SNP's callset.
gatk VariantFiltration --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $VFVCFD/selected_raw_snp_$run.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_SNP_filter" \
-O $VFVCFD/filtered_SNP_data_$run.vcf


### INDELs
echo "Extract and filter INDEL's from the call set."

	#3. Extract the INDELS from the ORIGINAL call set.
gatk SelectVariants --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $GVCFD/genotyped_data_$run.vcf \
--select-type-to-include INDEL \
-O $VFVCFD/selected_raw_indels_$run.vcf


	#4.Apply the filters to the INDEL's callset.
gatk VariantFiltration --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-V $VFVCFD/selected_raw_indels_$run.vcf \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "my_INDEL_filter" \
-O $VFVCFD/filtered_INDEL_data_$run.vcf


# 5. Combine Variants after using SNPS and INDELS filtering into a single file and get it ready for annotation.
gatk MergeVcfs --TMP_DIR=$TMP \
-R $HG19/ucsc.hg19.fasta \
-I $VFVCFD/filtered_SNP_data_$run.vcf \
-I $VFVCFD/filtered_INDEL_data_$run.vcf \
-O $VFVCFD/filtered_INDEL_SNP_data_$run.vcf

echo '\nHard filtering DONE' 










echo "\n\n\n- GENOTYPE POSTERIOR RECALCULATION BASED ON FAMILY PRIOR INFORMATION (GATK) "
echo "-------------------------------------------------------------------------------\n"


if [ "$ped" != "null" ]; then

	# Recalculate genotype posterior probabilities based on family prior information.
	# Compute the most likely genotype combination of trios and parent/child pairs given their genotype likelihoods and a mutation prior;
	# Phase all sites were parent/child transmission can be inferred unambiguously, ¿DO WE WANT THAT INFO?

	gatk CalculateGenotypePosteriors \
	   -R $HG19/ucsc.hg19.fasta \
	   -V $VFVCFD/filtered_INDEL_SNP_data_$run.vcf \
	   -ped $ped \
	   -O $VFVCFD/filtered_INDEL_SNP_data_GP_$run.vcf \
	   --skip-population-priors  

	# Filter variants based on GQ: genotypes with GQ < 20 based on the posteriors are flagged for posterior filtering. 
	gatk VariantFiltration \
	-R /mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle/ucsc.hg19.fasta \
	-V $VFVCFD/filtered_INDEL_SNP_data_GP_$run.vcf \
	--filter-expression "GQ < 20" \
	--filter-name "lowGQ" \
	-O $VFVCFD/filtered_INDEL_SNP_data_GP_GQfiltered_$run.vcf 


	# Annotation of High and Low Confidence De Novo mutations: 
	gatk3 --analysis_type VariantAnnotator \
	-R $HG19/ucsc.hg19.fasta \
	-V $VFVCFD/filtered_INDEL_SNP_data_GP_GQfiltered_$run.vcf  \
	-o $VFVCFD/filtered_INDEL_SNP_data_GP_GQfiltered_DeNovoAnnot_$run.vcf \
	-A PossibleDeNovo  \
	-ped $ped 


	VCF_IN="${VFVCFD}/filtered_INDEL_SNP_data_GP_GQfiltered_DeNovoAnnot_$run.vcf"

	echo '\nPosterior genotyping and de novo mutations annotation DONE\n'


else


	VCF_IN="${VFVCFD}/filtered_INDEL_SNP_data_${run}.vcf"


fi













echo "\n\n\n- VARIANT ANNOTATION (VEP - ENSEMBL) "
echo "----------------------------------------\n"



mkdir $VEPVCFAD
echo "mkdir vep_vcf_annotated_data"

VCF_OUT="${VEPVCFAD}/vep_qual_filters_${run}.vcf"
VCF_FILTERED="${VEPVCFAD}/vep_qual_filters_annotated_${run}.vcf"
VCF_FILTERED_2="${VEPVCFAD}/vep_qual_filters_annotated_can_conseq_${run}.vcf"
VCF_FILTERED_3="${VEPVCFAD}/vep_qual_filters_annotated_can_conseq_freqs_${run}.vcf"
VCF_FINAL="${VEPVCFAD}/FinalVariantAnnotation_${run}.txt"





echo '\nFiltering by quality and chromosome...\n'

perl $FILTER_VEP \
-i $VCF_IN -o $VCF_OUT \
--filter "QUAL > 100 and FILTER = PASS" \
--filter "CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chrX,chrY" --force_overwrite 




echo '\nVEP annotation...\n'

perl $VEP \
--cache --offline --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR --v --fork 16 --assembly GRCh37 --fasta $VEP_FASTA --force_overwrite \
--biotype --regulatory --protein --symbol --allele_number --numbers --domains --uniprot --variant_class \
--canonical --vcf \
--sift p --polyphen p --af --max_af \
--format vcf \
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
--filter "Consequence is 3_prime_UTR_variant or Consequence is 5_prime_UTR_variant or Consequence is intron_variant or Consequence is splice_donor_variant or Consequence is splice_acceptor_variant or Consequence is splice_region_variant or Consequence is synonymous_variant or Consequence is missense_variant or Consequence is inframe_deletion or Consequence is inframe_insertion or Consequence is stop_gained or Consequence is frameshift_variant or Consequence is coding_sequence_variant" --force_overwrite 

perl $FILTER_VEP \
-i $VCF_FILTERED_2 -o $VCF_FILTERED_3 \
--filter "AF < 0.01 or not AF" \
--filter "(ExAC_EAS_AF < 0.01 or not ExAC_EAS_AF) and (1000Gp3_AF < 0.01 or not 1000Gp3_AF) and (1000Gp3_EUR_AF < 0.01 or not 1000Gp3_EUR_AF ) and (gnomAD_exomes_AF < 0.01 or not gnomAD_exomes_AF) and (ExAC_AF < 0.01 or not ExAC_AF) and (ExAC_Adj_AF < 0.01 or not ExAC_Adj_AF) and (ExAC_NFE_AF < 0.01 or not ExAC_NFE_AF) and (gnomAD_exomes_NFE_AF < 0.01 or not gnomAD_exomes_NFE_AF) and (gnomAD_genomes_AF < 0.01 or not gnomAD_genomes_AF)" --force_overwrite 

echo '\nVEP annotation and filtering DONE\n'



# Remove files 

echo '\n'
rm -r $TMP





echo "\n\n\n- OUTPUT PROCESSING "
echo "-------------------------\n"


utilitiesPath=$(dirname "$0")


echo '\nFrom VEP annotated VCF to txt file...\n'


python $utilitiesPath/vep2tsvlorena.py $VCF_FILTERED_3 \
-o $VCF_FINAL -t -g -f GT 






# echo "\n\n\n- Organise final VEP file "
# echo "---------------------------------\n"


# utilitiesPath=dir=$(dirname "$0")





# if [ "$ped" != "null" ]; then


# 	gatk VariantsToTable \
# 	-R $HG19/ucsc.hg19.fasta \
# 	-V $VCF_FILTERED_3 \
# 	-F CHROM -F POS -F TYPE -F REF -F ALT -F HET -F HOM-VAR  -F HOM-REF -GF GT -GF GQ -F hiConfDeNovo -F lowConfDeNovo \
# 	--split-multi-allelic -O $MDAP/SNV_filtAnnot.txt --add-output-vcf-command-line

# else


# 	gatk VariantsToTable \
# 	-R $HG19/ucsc.hg19.fasta \
# 	-V $VCF_FILTERED_3 \
# 	-F CHROM -F POS -F TYPE -F REF -F ALT -F HET -F HOM-VAR  -F HOM-REF -GF GT -GF GQ -F hiConfDeNovo -F lowConfDeNovo \
# 	--split-multi-allelic -O $MDAP/SNV_filtAnnot.txt --add-output-vcf-command-line

# fi











# Trio phasing

#  This tool performs two functions:

#     Compute the most likely genotype combination of trios and parent/child pairs given their genotype likelihoods and a mutation prior;

# The tool ultimately reports the genotype combination (and hence phasing) probability.
# Ambiguous sites are:

#     Sites where all individuals are heterozygous
#     Sites where there is a Mendelian violation

# Missing genotypes are handled as follows:

#     In parent/child pairs: If an individual genotype is missing at one site, the other one is phased if it is homozygous. No phasing probability is emitted.
#     In trios: If the child is missing, parents are treated as separate individuals and phased if homozygous. No phasing probability is emitted.
#     In trios: If one of the parents is missing, it is handled like a parent/child pair. Phasing is done unless both the parent and child are heterozygous and a phasing probability is emitted.
#     In trios: If two individuals are missing, the remaining individual is phased if it is homozygous. No phasing probability is emitted.




