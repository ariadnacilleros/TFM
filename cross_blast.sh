#!/bin/bash

for i in `seq 1 2 41615`;
do
	j=$(($i+1))
	sed -n "$i","$j"p blast_seq.fasta > ilmn.fasta
	bin/ncbi-blast-2.9.0+/bin/blastn -query ilmn.fasta -taxids 9606 -db /home/ariadna/bin/blastdb/refseq_rna -out example.txt 
	count=`grep ">" example.txt | sed 's/.*Homo//' | sed 's/),.*$//' |  sort | uniq | wc -l `
	echo $count
	if [ "$count" != "1" ]; then
		grep "Query=" example.txt >> cross_hyb.txt
	fi
done
