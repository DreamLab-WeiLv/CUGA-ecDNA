${samtools} view -b -F 1294 ${tumor}.rmdup.sort.bam > ${tumor}.discordants.unsorted.bam
${samtools} sort -@ 12 -o ${tumor}.discordants.bam ${tumor}.discordants.unsorted.bam
${samtools} view -h ${tumor}.rmdup.sort.bam | ${LUMPY_DIR}/lumpy_scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${tumor}.splitters.unsorted.bam
${samtools} sort -@ 12 -o ${tumor}.splitters.bam ${tumor}.splitters.unsorted.bam
${LUMPY_DIR}/lumpyexpress -B ${tumor}.rmdup.sort.bam,${normal}.rmdup.sort.bam -S ${tumor}.splitters.bam,${normal}.splitters.bam -D ${tumor}.discordants.bam,${normal}.discordants.bam -o lumpy_${tumor}.vcf
${LUMPY_DIR}/svtyper -B ${tumor}.rmdup.sort.bam,${normal}.rmdup.sort.bam -l ${tumor}.bam.json -i lumpy_${tumor}.vcf > lumpy_${tumor}.sv.vcf
