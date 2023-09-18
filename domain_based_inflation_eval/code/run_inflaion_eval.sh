# we have to erase empty lines in tblout from pfam in order to rhmmer (R) read the files
sed -i '/^[[:space:]]*$/d' ../../pfam/results/*.tblout
# run Rscript
