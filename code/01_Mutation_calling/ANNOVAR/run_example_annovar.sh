## 01. Run ANNOVAR
perl ${annovar}/table_annovar.pl \\
${prx}.pass.snp.indel.vcf \\
${annovar}/humandb \\
--outfile ${outdir} \\
--buildver hg38 \\
--protocol refGene,ensGene,cytoBand,avsnp150,clinvar_20220320,cosmic70,dbnsfp42a,1000g2015aug_all,exac03,gnomad30_genome,gnomad_exome,icgc28 \\
--operation g,g,r,f,f,f,f,f,f,f,f,f \\
--vcfinput \\
--thread 6 \\
--dot2underline \\
--nastring . \\
--remove

## 02. Filter 
awk -F "\t" 'BEGIN{OFS="\t"} NR==1 {print $0}' ${prx}.hg38_multianno.txt > ${prx}.title.txt
awk -F "\t" 'BEGIN{OFS="\t"} $131<=0.01 && $132<=0.01 || $22=="Conflicting_interpretations_of_pathogenicity" || $22=="Pathogenic" || $22=="Likely_pathogenic" {print $0}' ${prx}.hg38_multianno.txt > ${prx}.filtered.txt
cat ${prx}.title.txt ${prx}.filtered.txt > ${prx}.filtered.txt
gzip ${prx}.filtered.txt
mv ${prx}.filtered.txt.gz ${prx}.hg38_multianno.txt.gz

## 03. VCF to MAF
perl maf.pl $prx_list $indir $outdir hg38_multianno.txt.gz
