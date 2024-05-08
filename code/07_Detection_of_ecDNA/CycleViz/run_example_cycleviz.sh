${AmpliconSuite-pipeline-1.0.0-pipeline}/scripts/CAMPER.py --graph ${graph_data} --runmode bulk
${CV_SRC}/CycleViz.py --ref GRCh38 --cycles_file ${cycles_data} --gene_subset_file "BUSHMAN" --cycle ${cycle} -g ${graph_data} --rotate_to_min --figure_size_style normal --outname ${prx}
