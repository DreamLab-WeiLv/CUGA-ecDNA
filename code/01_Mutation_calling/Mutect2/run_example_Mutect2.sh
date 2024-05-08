## 01. Running BQSR
${gatk} --java-options "-Xmx30g" BaseRecalibrator \
        -R hg38.fa \
        -I ${prx}.rmdup.sort.bam \
        -O ${prx}.recal_data.table \
        --known-sites ${dbsnp} \
        --known-sites ${dbsnp1000G} \
        --known-sites ${dbindel1000G}

${gatk} --java-options "-Xmx30g" ApplyBQSR \
        -R hg38.fa \
        -I ${prx}.rmdup.sort.bam \
        -O ${prx}.rmdup.sort.bqsr.bam \
        --bqsr-recal-file ${prx}.recal_data.table
## 02. Somatic calling
${gatk} --java-options "-Xmx30g" Mutect2 \
  -R hg38.fa \
  -I ${tumor_prx}.rmdup.sort.bqsr.bam \
  -tumor ${tumor_prx} \
  -I ${normal_prx}.rmdup.sort.bqsr.bam \
  -normal ${normal_prx} \
  --germline-resource somatic-hg38_af-only-gnomad.hg38.vcf.gz \
  -max-mnp-distance 0  \
  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY \\
  -O ${prx}.unfilter.vcf.gz

${gatk} --java-options "-Xmx30g" FilterMutectCalls \
  -R hg38.fa \
  -V ${prx}.unfilter.vcf.gz \
  -O ${prx}.filter.vcf.gz 

## 03. Split muti-allelic sites
${gatk} --java-options "-Xmx30g" LeftAlignAndTrimVariants \
 -V ${prx}.filter.vcf.gz \
 -R hg38.fa \
 -no-trim true --split-multi-allelics true \
 -O ${prx}.filter.m2.SplitMulti.vcf

## 04. Sprate SNP and INDEL
${gatk} --java-options "-Xmx30g" SelectVariants \
 -V ${prx}.filter.m2.SplitMulti.vcf \
 -R hg38.fa \
 --select-type-to-include SNP \
 -O ${prx}.m2.SplitMulti.snp.vcf

${gatk} --java-options "-Xmx30g" SelectVariants \
 -V ${prx}.filter.m2.SplitMulti.vcf \
 -R hg38.fa \
 --select-type-to-include INDEL \
 -O ${prx}.m2.SplitMulti.indel.vcf
