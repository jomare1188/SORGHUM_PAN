#PBS -q par16
#PBS -l nodes=1:ppn=16

cd $PBS_O_WORKDIR
conda activate TRANSDECODER

ls ./../../assemblies/ | grep .fasta > ./../../assemblies/transcriptomes_list

rm transdecoder_commands
for i in $(cat ./../../assemblies/transcriptomes_list)
do
	GENOTYPE=$(echo ${i} | cut -f2 -d"_")
	echo TransDecoder.LongOrfs -S -t ./../../assemblies/$i --output_dir ./../results/$GENOTYPE >> transdecoder_commands
done

cat transdecoder_commands | parallel -j16
