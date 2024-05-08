export MOSEKLM_LICENSE_FILE=${MOSEK_DIR}
export AA_DATA_REPO=${data_repo}
export AA_SRC=${AmpliconArchitect-1.3.r6}/src
export PATH=${MOSEK_DIR}/8/tools/platform/linux64x86/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MOSEK_DIR}/8/tools/platform/linux64x86/bin
export AC_SRC=${AmpliconClassifier-1.0.0}
python3 ${AmpliconSuite-pipeline-1.0.0}/PrepareAA.py -s ${tumor} -t 10 --ref GRCh38 --sorted_bam ${tumor_bam} --normal_bam ${normal_bam} --cnvkit_dir ${cnvkit_PATH} -o ${tumor} 
python3 ${AmpliconArchitect-1.3.r6}/src/AmpliconArchitect.py --bed ${tumor}/${tumor}_AA_CNV_SEEDS.bed --bam ${tumor_bam} --out ${tumor} --ref GRCh38

