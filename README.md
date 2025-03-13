# A Standardized Workflow for Pathogenic Variant Identification
This repository contains a production-grade bioinformatics workflow for clinical whole exome sequencing (WES) analysis.
The pipeline demonstrates the successful identification of a novel pathogenic variant in _SLC19A3_ causing thiamine transporter-2 deficiency.

# Table of Contents

1. Overview
2. System Requirements
3. Installation
4. Directory Structure
5. Pipeline Components
6. Usage
7. Example Analysis
8. Validation

# Overview
This pipeline implements the WES workflow used in identifying the novel c.871G>C (p.Gly291Arg) variant in the _SLC19A3_ gene, as described in "Early Infantile Thiamine Transporter-2 Deficiency with Epileptic Spasms—A Phenotypic Spectrum with a Novel Mutation" (Mishra et al., 2021), I am a co-author in this paper. The workflow enables sensitive detection of rare pathogenic variants in clinical samples, with particular emphasis on homozygous variants in consanguineous families.

# System Requirements
Hardware Requirements


## Exome Sequencing Pipeline

This pipeline processes whole exome sequencing data from raw reads to variant calling.
### 1. Set up and configuration

### 2. Quality Control Step

```bash
# Count reads in fastq files
zcat *R1*.fastq.gz | echo $((`wc -l`/4))
zcat *R2*.fastq.gz | echo $((`wc -l`/4))

# Run FastQC on all fastq files
while read line; do
    fastqc -t 8 "$line"	
done < fastq_list.txt
```



Counts reads in paired-end fastq files to ensure data completeness

Runs FastQC to assess sequence quality (base quality scores, GC content, sequence duplication, adapter content)

Generates HTML reports for quality visualization


### 3. Alignment to Reference Genome

```bash
# Align reads to reference genome
bwa mem -t 12 -v 1 /$HGDATA/hg19/hg19.fa -M -R "@RG\tID:Patient\tPL:illumina\tSM:SGRH" "$outdir"/*R1*.gz "$outdir"/*R2*.gz > "$outdir"/RAW_fastq_mem.sam
```
Uses BWA-MEM algorithm to align paired-end reads to hg19 reference genome

-t 12 uses 12 threads for faster processing

-M flags split hits as secondary (for compatibility with Picard)

-R adds read group information (required for GATK)

Output is in SAM format


### 4. Post-Alignment Processing

```bash
# Filter multimapped reads
awk '{if ($7!="*") print $0}' RAW_fastq_mem.sam > RAW_fastq_mem_multimapped_unmapped_filtered.sam

# Sort SAM file
samtools sort -m 8G -@ 4 -o RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup.sam RAW_fastq_mem_multimapped_unmapped_filtered.sam

# Convert to BAM format and index
samtools view -bhS RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup.sam > RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup.bam
samtools index RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup.bam

# Clean up SAM files to save disk space
rm -f *.sam
```
Filters out multimapped reads (improves specificity)

Sorts reads by genomic coordinates

Converts SAM to compressed BAM format

Creates index for random access to BAM file

### 5. GATK Pre-processing for Improved Accuracy
```bash
# Identify targets for realignment
java -Xmx8g -jar /$PACKAGE/GATK/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /$HGDATA/hg19/hg19.fa -o intervalsList.intervals -I RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup.bam -known /$RESOURCES/variant_calling_data/1000g_gold_standard.indels.hg19.sites.vcf

# Perform local realignment around indels
java -Xmx8g -jar /$PACKAGE/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R /$HGDATA/hg19/hg19.fa -I RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup.bam -targetIntervals intervalsList.intervals -known /$RESOURCES/variant_calling_data/1000g_gold_standard.indels.hg19.sites.vcf -o RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned.bam --filter_bases_not_stored
```
RealignerTargetCreator identifies regions needing realignment (usually around indels)

IndelRealigner performs local realignment to minimize mismatches around indels

Reduces false positive variant calls, especially around indels

### 6. Base Quality Score Recalibration (BQSR)
```bash
# Generate recalibration table
java -Xmx8g -jar /$PACKAGE/GATK/GenomeAnalysisTK.jar -I RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned.bam -R /$HGDATA/hg19/hg19.fa -T BaseRecalibrator -knownSites /$RESOURCES/variant_calling_data/dbsnp_sorted.hg19.vcf -o RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned.bam.pre.recal.table

# Apply recalibration
java -Xmx8g -jar /$PACKAGE/GATK//GenomeAnalysisTK.jar -T PrintReads -R /$HGDATA/hg19/hg19.fa -I RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned.bam -BQSR RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned.bam.pre.recal.table -o RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned_recal.bam
```
BaseRecalibrator adjusts quality scores using machine learning based on known variants

PrintReads applies these adjustments to create a new BAM file

Improves accuracy of variant calling by correcting systematic errors in base quality scores

### 7. Variant Calling with HaplotypeCaller
```bash
# Call variants using HaplotypeCaller
java -Xmx8g -jar /$PACKAGE/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /$HGDATA/hg19/hg19.fa -I RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned_recal.bam -o RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned_recal.bam.HC.vcf --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30
```
Uses GATK's HaplotypeCaller to identify SNPs and indels

Creates local haplotypes to improve variant calling accuracy

Sets confidence thresholds for variant emission (10) and calling (30)

Outputs a VCF file with raw variant calls

### 8. Variant Quality Score Recalibration (VQSR) - SNPs
```bash
# Build recalibration model for SNPs
java -Xmx4g -jar /$PACKAGE/GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R /$HGDATA/hg19/hg19.fa -input RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned_recal.bam.HC.vcf -resource:hapmap,VCF,Known=false,training=true,truth=true,prior=15.0 /$RESOURCES/variant_calling_data/hapmap_3.3.hg19.sites_sorted.vcf -resource:omni,VCF,Known=false,training=true,truth=false,prior=12.0 /$RESOURCES/variant_calling_data/1000G_omni2.5.hg19.sites_sorted.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=6.0 /$RESOURCES/variant_calling_data/dbsnp_sorted.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 /$RESOURCES/variant_calling_data/1000G_phase1.snps.high_confidence.hg19.vcf -an QD -an DP -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.tranches.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R

# Apply recalibration to SNPs
java -Xmx8g -jar /$PACKAGE/GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R /$HGDATA/hg19/hg19.fa --input RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned_recal.bam.HC.vcf --mode SNP --ts_filter_level 99.0 -recalFile recalibrate_SNP.tranches.recal -tranchesFile recalibrate_SNP.tranches -o RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned_recal_snps_raw_indels.vcf
```
VariantRecalibrator builds a Gaussian mixture model based on known variants (HapMap, 1000G)

Uses multiple annotation metrics (QD, DP, MQ, etc.) to assess variant quality

ApplyRecalibration applies this model to filter SNPs with 99% sensitivity threshold

Separates high-quality variants from potential false positives

### 9. Variant Quality Score Recalibration (VQSR) - Indels
```bash
# Build recalibration model for indels
java -Xmx8g -jar /$PACKAGE/GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R /$HGDATA/hg19/hg19.fa -input RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned_recal_snps_raw_indels.vcf -resource:1000G,VCF,Known=false,training=true,truth=true,prior=10.0 /$RESOURCES/variant_calling_data/1000g_gold_standard.indels.hg19.sites.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=6.0 /$RESOURCES/variant_calling_data/dbsnp_sorted.hg19.vcf -an QD -an DP -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_INDEL.tranches.recal -tranchesFile recalibrate_INDEL.tranches --maxGaussians 4 -rscriptFile recalibrate_INDEL_plots.R

# Apply recalibration to indels
java -Xmx8g -jar /$PACKAGE/GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R /$HGDATA/hg19/hg19.fa --input RAW_fastq_mem_multimapped_unmapped_filtered_sorted_rmdup_realigned_recal_snps_raw_indels.vcf --mode INDEL --ts_filter_level 99.0 -recalFile recalibrate_INDEL.tranches.recal -tranchesFile recalibrate_INDEL.tranches -o recalibrated_SNP_indels.vcf
```
Similar to SNP recalibration but tailored for indels

Uses fewer Gaussians (4) since indel metrics have different distributions

Produces a final VCF with high-quality SNPs and indels
