base_dir=$PWD
orthogroups_dir=/home/dmpachon/sorghum_orthofinder/results/I2/Results_Oct24/Orthogroup_Sequences
cd $orthogroups_dir
rm -f ${base_dir}/sequences_per_orthogroup.csv
count=0
for i in $(ls)
do
	lines=$(grep -c ">" ${i})
	echo ${i},${lines} >> ${base_dir}/../results/sequences_per_orthogroup.csv
	((count++))
	echo $count
done
cd ${base_dir}

awk -F ',' '{if ($2 > 4) print $0}' ${base_dir}/../results/sequences_per_orthogroup.csv | cut -f1 -d, > ${base_dir}/../results/sequences_per_orthogroup_filtered4444.csv

sed -i ${base_dir}/../results/sequences_per_orthogroup_filtered4.csv 's/.fa//g'
