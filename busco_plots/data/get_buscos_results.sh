for i in $(cat genotypes)
do
	cp -v ../../busco/results/${i}_busco/short_summary.specific.embryophyta_odb10.${i}_busco.json ${i}_busco.json
done


