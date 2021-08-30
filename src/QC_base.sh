basedir=$1 # gwas summary statistics of base populations
basepop=$2

#mkdir ${cur}/output
#mkdir ${cur}/output/base
#output=${cur}/output/base
awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' ${basedir}/${basepop}.gwas.txt|\
	gzip  > ${basedir}/${basepop}.gz

gunzip -c ${basedir}/${basepop}.gz |\
	awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
	gzip - > ${basedir}/${basepop}.nodup.gz

gunzip -c ${basedir}/${basepop}.nodup.gz |\
	awk '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' |\
	gzip > ${basedir}/${basepop}.QC.gz

rm ${basedir}/${basepop}.gz ${basedir}/${basepop}.nodup.gz
