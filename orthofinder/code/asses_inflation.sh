## activate conda envirormment with mcl installed 
conda activate ORTHOFINDER
### MAKE TABLE WITH CLM-info
## you must have a file with the inflations values (each value one row and one decimal point)
#ls ./../results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder/ > ./../data/inflations.txt

# make a variable with the path to easier the things
results="./../results/results_2/Results_Jun14/WorkingDirectory/OrthoFinder"

# make directory results
mkdir -p ./../results/inflation

# remove files created in loops to avoid overwrite them
rm -f ./../results/inflation/clm-info_table.csv

# make the for loop to run clm info for each inflation value 
for i in $(cat ./../data/inflations.txt)
do
	echo $i
	clm info $results/${i}/WorkingDirectory/OrthoFinder_graph.txt  $results/${i}/WorkingDirectory/clusters_OrthoFinder_I${i}.txt | awk '{ for (i=1; i<=NF; i++) { split($i, arr, "="); printf "%s%s%s", sep, arr[2], (i==NF) ? "\n" : "," } }' >> ./../results/inflation/clm-info_table.csv
done

# put inflation value to the table
paste -d"," ./../data/inflations.txt ./../results/inflation/clm-info_table.csv > ./../results/inflation/clm-info_table2.csv

## now we make pretty table with clm info output for inflation values
# make the header of the table
echo "inflation,efficiency,massfrac,areafrac,source,clusters,max,ctr,avg,min,DGI,TWI,TWL,sgl,qrt" > ./../results/inflation/header_clm-info.csv
# put togheter header and table
cat ./../results/inflation/header_clm-info.csv ./../results/inflation/clm-info_table2.csv > ./../results/inflation/table_clm.csv

