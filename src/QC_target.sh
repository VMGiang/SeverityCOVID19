basedir=$1
basepop=$2
tardir=$3
tarpop=$4

#mkdir ${cur}/output/target
#tardir=${cur}/output/target

plink --bfile ${tardir}/${tarpop} \
	--maf 0.01 \
	--hwe 1e-6 \
    	--geno 0.01 \
	--mind 0.1 \
	--write-snplist \
	--make-just-fam \
	--out ${tardir}/${tarpop}.QC

plink --bfile ${tardir}/${tarpop} \
	--keep ${tardir}/${tarpop}.QC.fam \
	--extract ${tardir}/${tarpop}.QC.snplist \
	--indep-pairwise 200 50 0.25 \
	--out ${tardir}/${tarpop}.QC

plink --bfile ${tardir}/${tarpop} \
	--extract ${tardir}/${tarpop}.QC.prune.in \
	--keep ${tardir}/${tarpop}.QC.fam \
	--het \
	--out ${tardir}/${tarpop}.QC

Rscript --no-save calculation_F_score.R ${tardir} ${tarpop}

Rscript --no-save filter_mismatching_snps.R ${tardir} ${basedir} ${tarpop} ${basepop}

plink 	--bfile ${tardir}/${tarpop} \
	--extract ${tardir}/${tarpop}.QC.prune.in \
	--keep ${tardir}/${tarpop}.valid.sample \
	--check-sex \
	--out ${tardir}/${tarpop}.QC
Rscript --no-save check_sex.R ${tardir} ${tarpop} 

plink	--bfile ${tardir}/${tarpop} \
	--extract ${tardir}/${tarpop}.QC.prune.in \
	--keep ${tardir}/${tarpop}.QC.valid \
	--rel-cutoff 0.125 \
	--out ${tardir}/${tarpop}.QC
plink \
	--bfile ${tardir}/${tarpop} \
	--make-bed \
	--keep ${tardir}/${tarpop}.QC.rel.id \
	--out ${tardir}/${tarpop}.QC \
	--extract ${tardir}/${tarpop}.QC.snplist \
	--exclude ${tardir}/${tarpop}.mismatch \
	--a1-allele ${tardir}/${tarpop}.a1






