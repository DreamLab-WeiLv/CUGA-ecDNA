 ${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py --normalBam ${normal}.rmdup.sort.bam --tumorBam ${tumor}.rmdup.sort.bam --referenceFasta hg38.fa --indelCandidates ${tumor}_candidateSmallIndels.vcf.gz --runDir ${tumor}
