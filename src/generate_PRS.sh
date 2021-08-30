basedir=$1
basepop=$2
tardir=$3
tarpop=$4

Rscript Transform_base.R ${basedir} ${basepop}

plink --bfile ${tardir}/${tarpop}.QC \
	--clump-p1 1 \
	--clump-r2 0.1 \
	--clump-kb 250 \
	--clump ${basedir}/${basepop}.QC.Transformed \
	--clump-snp-field SNP \
	--clump-field P \
	--out ${tardir}/${tarpop}

awk 'NR!=1{print $3}' ${tardir}/${tarpop}.clumped >  ${tardir}/${tarpop}.valid.snp

awk '{print $3,$8}' ${basedir}/${basepop}.QC.Transformed > ${tardir}/SNP.pvalue

touch ${tardir}/range_list
for range_list in "0.001 0 0.001" "0.05 0 0.05" "0.1 0 0.1" "0.2 0 0.2" "0.3 0 0.3" "0.4 0 0.4" "0.5 0 0.5"
do
  echo $range_list >> ${tardir}/range_list
done


mkdir ${tardir}/PRS

plink \
    --bfile ${tardir}/${tarpop}.QC \
    --score ${basedir}/${basepop}.QC.Transformed 3 4 12 header \
    --q-score-range ${tardir}/range_list ${tardir}/SNP.pvalue \
    --extract ${tardir}/${tarpop}.valid.snp \
    --out ${tardir}/PRS/${tarpop}

#rm ${tardir}/*.log 
#rm ${tardir}/range_list ${tardir}/*.a1 ${tardir}/*.mismatch ${tardir}/*.valid.sample ${tardir}/*.QC.*




