## 01. Merge Mutect2 and strelka2 Indel vcf 
${bedtools} intersect -u -a ${Mutect2_Vcf} -b ${strelka_Vcf} > ${prx}_merge.indel.vcf

## 02. Merge Indel and SNP vcf 
cat ${Mutect2}.snp.vcf ${prx}_merge.indel.vcf > ${prx}.snp.indel.vcf

## 03. Filter PASS
${gatk} --java-options "-Xmx40g" SelectVariants -R ${ref} -V ${prx}.snp.indel.vcf --exclude-filtered -O ${prx}.pass.snp.indel.vcf
