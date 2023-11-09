## activate conda envirormment with mcl installed 
module load mcl/14-137

### MAKE TABLE WITH CLM-info
# you must have a file with the inflations values (each value one row and one decimal point) in ../data/inflations.txt

# remove files created in loops to avoid overwrite them
rm -f ../results/clm-info_table.csv commands

# make the for loop to run clm info for each inflation value 
for i in $(cat ./../data/inflations.txt)
do
	echo $i
	clm info ../data/graphs/OrthoFinder_graph_I${i}.txt  ../data/clusters/clusters_OrthoFinder_I${i}.txt | awk '{ for (i=1; i<=NF; i++) { split($i, arr, "="); printf "%s%s%s", sep, arr[2], (i==NF) ? "\n" : "," } }' >> ../results/clm-info_table.csv
done

# put inflation value to the table
paste -d"," ./../data/inflations.txt ./../results/clm-info_table.csv > ./../results/clm-info_table2.csv

## now we make pretty table with clm info output for inflation values
# make the header of the table
echo "inflation,efficiency,massfrac,areafrac,source,clusters,max,ctr,avg,min,DGI,TWI,TWL,sgl,qrt" > ./../data/header_clm-info.csv
# put togheter header and table
cat ./../data/header_clm-info.csv ./../results/clm-info_table2.csv > ./../results/table_clm.csv

rm -f ./../results/clm-info_table2.csv ./../results/clm-info_table.csv ./../data/header_clm-info.csv

## clm dist
#clm dist --progress -o ./../results/clm-dist.tsv ../data/clusters/clusters_OrthoFinder_I*
## clm meet
#clm meet -o ./../results/clm-meet.tsv ../data/clusters/clusters_OrthoFinder_I*
# clm meet 2?
#clm dist ./../results/clm-meet.tsv ../data/clusters/clusters_OrthoFinder_I* > ./../results/clm-dist-with_meet.tsv

#clm dist --progress -o ./../results/inflation/clm-dist.tsv $results/1.1/WorkingDirectory/clusters_OrthoFinder_I1.1.txt $results/1.3/WorkingDirectory/clusters_OrthoFinder_I1.3.txt $results/1.5/WorkingDirectory/clusters_OrthoFinder_I1.5.txt $results/1.8/WorkingDirectory/clusters_OrthoFinder_I1.8.txt $results/2.0/WorkingDirectory/clusters_OrthoFinder_I2.0.txt $results/2.5/WorkingDirectory/clusters_OrthoFinder_I2.5.txt $results/3.0/WorkingDirectory/clusters_OrthoFinder_I3.0.txt $results/4.0/WorkingDirectory/clusters_OrthoFinder_I4.0.txt $results/5.0/WorkingDirectory/clusters_OrthoFinder_I5.0.txt $results/6.0/WorkingDirectory/clusters_OrthoFinder_I6.0.txt $results/7.0/WorkingDirectory/clusters_OrthoFinder_I7.0.txt $results/10.0/WorkingDirectory/clusters_OrthoFinder_I10.0.txt $results/15.0/WorkingDirectory/clusters_OrthoFinder_I15.0.txt $results/20.0/WorkingDirectory/clusters_OrthoFinder_I20.0.txt

##
