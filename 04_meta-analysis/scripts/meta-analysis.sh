
infiles=$1
outfile=$2
metal_in=$3

# infiles="/newhome/tb13101/PACE/pace_ad/data/meta_analysis/m1a/ALSPAC_AD_m1a.txt /newhome/tb13101/PACE/pace_ad/data/meta_analysis/m1a/ad_childhood_sex_sv_NoCell_model1a_GOYA_20210430.txt.gz"
# outfile="results/metal-res/m1a.txt"
# metal_in="/newhome/tb13101/PACE/pace_ad/meta-analysis/data/metal-infiles/m1a.txt"

metal_outfile=${outfile%.txt}
metal="/newhome/tb13101/PACE/random-metal/executables/metal"


IFS=' ' read -r -a infiles_array <<< "$infiles"

echo ${infiles_array[0]}

rm -f ${metal_in}
touch ${metal_in}

echo "OUTFILE ${metal_outfile} .txt" >> ${metal_in}
echo "" >> ${metal_in}
echo "MARKER probeID" >> ${metal_in}
echo "EFFECTLABEL BETA" >> ${metal_in}
echo "STDERRLABEL SE" >> ${metal_in}
echo "PVALUELABEL P_VAL" >> ${metal_in}
echo "SEPARATOR TAB" >> ${metal_in}
echo "COLUMNCOUNTING LENIENT" >> ${metal_in}
echo "MINMAXFREQ OFF" >> ${metal_in}
echo "AVERAGEFREQ OFF" >> ${metal_in}
echo "GENOMICCONTROL OFF" >> ${metal_in}
echo "USESTRAND OFF" >> ${metal_in}
echo "SCHEME STDERR" >> ${metal_in}
echo "" >> ${metal_in}


for file in "${infiles_array[@]}"
do

	echo $file
	echo "PROCESS ${file}" >> ${metal_in}

done

echo "" >> ${metal_in}
echo "ANALYZE HETEROGENEITY" >> ${metal_in}
echo "CLEAR" >> ${metal_in}
echo "QUIT" >> ${metal_in}

${metal} ${metal_in}
mv ${metal_outfile}1.txt ${outfile}
