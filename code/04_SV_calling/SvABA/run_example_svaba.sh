${svaba} run -t ${tumor_bam} -n ${normal_bam} -a ${tumor} -p 10 -D dbsnp_indel.vcf -G hg38.fa
python svaba_transfer.py ${tumor}.svaba.somatic.sv.vcf > ${tumor}.svaba.sv.vcf
