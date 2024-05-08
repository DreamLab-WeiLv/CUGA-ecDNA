sh ${AC_DIR}/make_input.sh ${prx} ${prx}
export AA_DATA_REPO=${data_repo_DIR}
export AC_SRC=${AmpliconClassifier-1.0.0_DIR}
python ${AC_SRC}/amplicon_classifier.py --ref GRCh38 --input ${input} --report_complexity --verbose_classification -o ${prx}
