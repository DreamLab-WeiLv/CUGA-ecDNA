${RSEM_DIR}/rsem-calculate-expression --forward-prob 0.5 \
--paired-end \
-p 12 --no-bam-output \
--bam ${alignDir}/${prx}Aligned.toTranscriptome.out.bam \
${RSEM_DIR}/ref/human_gencode \
${OUT_DIR}
