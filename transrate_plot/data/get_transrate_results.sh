for i in $(cat genotypes)
do
	cp -v ../../transrate/results/${i}_transrate/assemblies.csv ${i}_tr.csv
done
