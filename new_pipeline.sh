#!/bin/bash
source /home/norrissw/bin/source_BNFO620.sh
source /usr/global/jdk7/java_env.sh
REF=/gpfs_fs/home/bnfo620/laahirie/hg19.reorder.fa
THREAD=15
F1=$1
F2=$2
VCFdb=/gpfs_fs/bnfo620/exome_data/dbsnp_138.hg19.vcf
VCFI=/gpfs_fs/bnfo620/exome_data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
VCF=/gpfs_fs/bnfo620/exome_data/hapmap_3.3.hg19.sites.vcf
#echo "Executing bwa mem."
#bwa mem -t $THREAD $REF $F1 $F2 > ${F1}.merged.sam
echo "Finished writing SAM."
#echo "Executing readgroups."
#java -jar /usr/global/blp/picard-tools-1.95/AddOrReplaceReadGroups.jar INPUT=${F1}.merged.sam OUTPUT=${F1}.RG.merged.sam RGLB=lib RGPL=illumina RGPU=flowcell-barcode.lane RGSM=sample
echo "Finished writing ReadGroup."
#echo "Executing samtools view."
#/usr/global/blp/bin/samtools view -bS ${F1}.RG.merged.sam -o ${F1}.merged.bam
echo "Finished writing BAM."
#echo "Executing picard-tools ReorderSam."
#java -jar /usr/global/blp/picard-tools-1.95/ReorderSam.jar I=${F1}.merged.bam o=${F1}.reorder.merged.bam R=$REF
echo "Finished reordering BAM."
#echo "Executing picart-tools SortSam."
#java -jar /usr/global/blp/picard-tools-1.95/SortSam.jar INPUT=${F1}.reorder.merged.bam OUTPUT=${F1}.sorted.merged.bam SORT_ORDER=coordinate
echo "Finished sorting reordered BAM."
#echo "Executing picard-tools MarkDuplicates."
#java -jar /usr/global/blp/picard-tools-1.95/MarkDuplicates.jar INPUT=${F1}.sorted.merged.bam OUTPUT=${F1}.dup_removed.merged.bam METRICS_FILE=${F1}.merged.metrics.txt
echo "Finished removing duplicates from BAM."
#echo "Executing picard-tools BuildBamIndex."
#java -jar /usr/global/blp/picard-tools-1.95/BuildBamIndex.jar INPUT=${F1}.dup_removed.merged.bam
echo "Finished indexing BAM."
#echo "Executing picard-tools BamIndexStats."
#java -jar /usr/global/blp/picard-tools-1.95/BamIndexStats.jar INPUT=${F1}.dup_removed.merged.bam > ${F1}.merged.bamStats.txt
echo "Finished generating statistics for BAM."
#echo "Executing GATK RealignerTargetCreator."
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I ${F1}.dup_removed.merged.bam -o ${F1}.merged.bam.target_creater.intervals
echo "Finished step 1 of realignment."
#echo "Executing GATK IndelRealigner."
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -I ${F1}.dup_removed.merged.bam -R $REF -T IndelRealigner -targetIntervals ${F1}.merged.bam.target_creater.intervals -o ${F1}.dup_removed.realigned.merged.bam
echo "Finished step 2 of realignment."
#echo "Executing picard-tools FixMateInformation."
#java -jar /usr/global/blp/picard-tools-1.95/FixMateInformation.jar INPUT=${F1}.dup_removed.realigned.merged.bam OUTPUT=${F1}.dup_removed.new.realigned.merged.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
echo "Finished step 3 of realignment."
#echo "Executing GATK BaseRecalibrator."
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I ${F1}.dup_removed.new.realigned.merged.bam -knownSites $VCFdb  -knownSites $VCFI -o ${F1}.merged.baseRe.grp
echo "Finished recalibrating BAM."
#echo "Executing GATK PrintReads."
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -T PrintReads -R $REF -I ${F1}.dup_removed.new.realigned.merged.bam -BQSR ${F1}.merged.baseRe.grp -o ${F1}.recal.merged.bam
echo "Finished filtering BAM."
#echo "Executing GATK HaplotypeCaller."
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF -I ${F1}.recal.merged.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o ${F1}.raw.variants.merged.g.vcf
echo "Finished calling SNPs and indels."
#echo "Executing GATK GenotypeGVCFs."
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -R $REF -T GenotypeGVCFs --variant ${F1}.raw.variants.merged.g.vcf -o output.merged.vcf
echo "Finished aggregating GVCFs."
#echo "Executing GATK VariantAnnotator."
#java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -R $REF -T VariantAnnotator -I ${F1}.recal.merged.bam -o output.fixed.merged.vcf -A Coverage --variant output.merged.vcf --dbsnp $VCFdb
echo "Finished annotating VCF."
echo "Executing GATK VariantRecalibrator."
java -jar /usr/global/blp/GenomeAnalysisTK-3.1.1/GenomeAnalysisTK.jar -R $REF -T VariantRecalibrator -input output.merged.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $VCF -resource:omni,known=false,training=true,truth=false,prior=12.0 $VCFI -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $VCFdb -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode SNP -recalFile output.cf.recal -tranchesFile output.cf.tranches -rscriptFile output.cf.plots.R
echo "Finished recalibrating VCF."






