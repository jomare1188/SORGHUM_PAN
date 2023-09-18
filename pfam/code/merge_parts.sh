cd /home/dmpachon/pfam_inflation_inflation/results/parts
for i in $(cat prefixes.list)
do
	ls | grep $i | grep pfam.tblout$ | xargs cat | grep -v "#" > $i.full_pfam.tblout.tmp
	ls | grep $i | grep pfam.out$ | xargs cat | grep -v "#" > ./../$i.full_pfam.out
#ls | grep $i | grep domtblout$ | xargs cat | grep -v "#" > $i.full_pfam.domtblout.tmp
	cat ./../../data/header_pfam.tblout $i.full_pfam.tblout.tmp > ./../$i.full_pfam.tblout && rm -f $i.full_pfam.tblout.tmp
#	cat ./../../data/header_pfam.out $i.full_pfam.out.tmp > ./../$i.full_pfam.out && rm -f $i.full_pfam.out.tmp
#cat ./../../data/header_pfam.domtblout $i.full_pfam.domtblout.tmp > ./../$i.full_pfam.domtblout && rm -f $i.full_pfam.domtblout.tmp
done
cd ./../../code/

