${delly} call -o ${tumor}.sv.bcf -g hg38.fa ${tumor}.rmdup.sort.bam ${normal}.rmdup.sort.bam
echo -e "${tumor}\ttumor\n${normal}\tcontrol" > ${tumor}.tsv
${delly} filter -f somatic -o ${tumor}.sv.ft.bcf -s ${tumor}.tsv ${tumor}.sv.bcf
${bcftools} view ${tumor}.sv.ft.bcf > ${tumor}.sv.vcf
awk -F "\t" 'BEGIN{OFS="\t"} $7 == "PASS" {print $0}' ${tumor}.sv.vcf > ${tumor}.delly.sv.vcf

