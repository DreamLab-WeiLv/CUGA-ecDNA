snakemake -s  Snakefile pretrained_predict -p -w 1000 --cores 5 --config mutations=${maf} output_dir=${prx_outdir} trained_classifier="2020plus_10k.Rdata"
