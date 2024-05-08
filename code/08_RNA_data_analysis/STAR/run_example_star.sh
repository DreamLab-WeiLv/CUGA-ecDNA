${STAR}/bin/Linux_x86_64_static/STAR \
--readFilesCommand zcat \
--quantMode TranscriptomeSAM GeneCounts \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped None \
--genomeDir ${STAR_DIR}/STAR_static \
--readFilesIn ${cleanDir}/r1_val_1.fq.gz r2_val_2.fq.gz \
--outFileNamePrefix ${alignDir}/${prx}
