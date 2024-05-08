${MANTA_DIR}/bin/configManta.py --normalBam ${normal}.rmdup.sort.bam --tumorBam ${tumor}.rmdup.sort.bam --referenceFasta hg38.fa --runDir ${tumor}
${currentdir}/${tumor}/runWorkflow.py -j 10
gzip -d -c ${tumor}/results/variants/somaticSV.vcf.gz > ${tumor}/results/variants/somaticSV.vcf
${MANTA_DIR}/libexec/convertInversion.py ${samtools} hg38.fa ${tumor}/results/variants/somaticSV.vcf > ${tumor}/results/variants/manta_${tumor}.sv.vcf
awk -F "\t" 'BEGIN{OFS="\t"} $7 == "PASS" {print $0}' ${tumor}/results/variants/manta_${tumor}.sv.vcf > ${tumor}.manta.sv.pass.vcf
