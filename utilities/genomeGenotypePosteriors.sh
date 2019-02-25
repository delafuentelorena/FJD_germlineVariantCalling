#!/bin/sh
export SFT=/mnt/genetica3/marius/pipeline_practicas_marius/software
alias bwa='$SFT/bwa/bwa'
alias samtools='$SFT/samtools/samtools'
alias picard='java -jar $SFT/picard/build/libs/picard.jar'
alias gatk='java -jar $SFT/gatk/build/libs/gatk-package-4.0.6.0-22-g9d9484f-SNAPSHOT-local.jar'
alias gatk3='java -jar /mnt/genetica3/GeneticaPipeDB_updated/gatk-3.8/GenomeAnalysisTK.jar'


# arguments

inputvep="/mnt/genetica3/lorena/TRIO/results/variant_filtration_vcf_data/filtered_INDEL_SNP_data_trio_2019_02_13_17_04_56.vcf"
output="trio_before_vep/"



# Recalculate genotype posterior probabilities based on family prior information.
gatk CalculateGenotypePosteriors \
	-R /mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle/ucsc.hg19.fasta \
   -V $inputvep \
   -ped /mnt/genetica2/NGS_data/Exomas-WES/CNAG-CONSYN-2013/Kabuki1/RESULTS-OLD-PIPELINE/kabuki1.ped \
   -O $output/output.withPosteriors.vcf \
   --skip-population-priors  ####  ¿?¿?



# Filter variants based on GQ: genotypes with GQ < 20 based on the posteriors are filtered out. GQ20 is widely accepted as a good threshold for genotype accuracy, indicating that there is a 99% chance that the genotype in question is correct.
gatk VariantFiltration \
-R /mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle/ucsc.hg19.fasta \
-V $output/output.withPosteriors.vcf \
--filter-expression "GQ < 20" \
--filter-name "lowGQ" \
-O $output/output.withPosteriors_filtered.vcf 
#Do we need to filter tagged genotypes?



# Annotation of High and Low Confidence De Novo mutations: 
gatk3 --analysis_type VariantAnnotator \
-R /mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle/ucsc.hg19.fasta \
-V $output/output.withPosteriors_filtered.vcf \
-o $output/output.withPosteriors_filtered.PossibleDeNovo.vcf \
-A PossibleDeNovo  \
-ped /mnt/genetica2/NGS_data/Exomas-WES/CNAG-CONSYN-2013/Kabuki1/RESULTS-OLD-PIPELINE/kabuki1.ped 
# Extract de high and low condident de novo mutations -> hiConfDeNovo=13-0262 (vector of childs)



# Summary table:
# table - option 1 - gatk
gatk VariantsToTable \
-R /mnt/genetica3/marius/pipeline_practicas_marius/hg19bundle/ucsc.hg19.fasta \
-V $output/output.withPosteriors_filtered.PossibleDeNovo.vcf \
-F CHROM -F POS -F TYPE -F REF -F ALT -F HET -F HOM-VAR  -F HOM-REF -GF GT -GF GQ -F hiConfDeNovo -F lowConfDeNovo -F CSQ \
-O $output/test_table2.txt

#  \-F lowGQ


# table - option 2 - plink
## Biallelic only variants ___ (removes null alternates hated by plink)
GenomeAnalysisTK  SelectVariants \
-R ${refseq} -V ${vcf} \
-restrictAllelesTo BIALLELIC \
-o bi${vcf}

## Convert to plink binary
plink --allow-extra-chr --allow-no-sex \
--vcf bi.${vcf}  \
--make-bed \
--out bi.${vcf}

## Convert plink binary to basic dosage
plink --allow-extra-chr --allow-no-sex \
--bfile bi.${vcf} \
--recode A-transpose \
--out bi.${vcf}.recode


#/share/apps/bio/gatk-4.0.11.0/gatk --java-options "-Xmx10g" SelectVariants -R /shared/resources/hgRef/hg38/Homo_sapiens_assembly38.fasta \
#-V /home/manolis/GATK4/3.WES_Illumina/germSNV/4.VCF/storage/raw_filtered_by_VQSR_by_on_target_regions.vcf -O prova.vcf -exclude-non-variants -sn "12-928" -sn "12-929" -sn "12-931"


