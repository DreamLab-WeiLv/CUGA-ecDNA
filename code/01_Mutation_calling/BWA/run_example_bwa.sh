${bwa} -t 10 -K 100000000 -R "@RG\tID:${sample_id}\tSM:${sample_id}\tLB:lib${sample_id}\tPL:illumina" -Y ${genome_ref} ${r1} ${r2} | ${samtools} view -@ 10 -1 - > ${prx}.bam
${gatk} --java-options MarkDuplicates --INPUT ${prx}.bam --OUTPUT ${prx}.rmdup.bam --METRICS_FILE ${prx}.metrics --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER "queryname"
${samtools} sort -@ 10 -o ${prx}.rmdup.sort.bam ${prx}.rmdup.bam
${samtools} index -@ 10 ${prx}.rmdup.sort.bam
