#!/bin/bash

## conda env create
mamba create -n anno_function -c bioconda -c conda-forge diamond==2.2.8 hmmer seqkit taxonkit csvtk pigz kobas kofamscan jcvi venn pixi
mamba activate anno_function
pip install pyhmmsearch

## NR (Using Magnoliopsida here, not all plants)
#prot.accession2taxid.gz: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
#taxdump.tar.gz: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
#NR: https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
taxonkit list --ids 3398 --indent "" > Magnoliopsida.taxid.txt
zcat prot.accession2taxid.gz | csvtk -t grep -f taxid -P Magnoliopsida.taxid.txt | csvtk -t cut -f accession.version > Magnoliopsida.taxid.acc.txt
seqkit grep -f Magnoliopsida.taxid.acc.txt nr > NR_Magnoliopsida.fa
diamond makedb --in NR_Magnoliopsida.fa --db NR_Magnoliopsida
grep -P '^>' NR_Magnoliopsida.fa | awk '/^>/{s=$0;sub(/^>/,"",s);match(s,/^[^ \t]+/);id=substr(s,RSTART,RLENGTH);desc=substr(s,RLENGTH+2);if(match(desc,/\[[^][]+\]/,m))print id"\t"substr(desc,1,RSTART+RLENGTH-1)}' > NR_Magnoliopsida_id2anno.tsv

## Uniprot
#uniprot_sport_plants.xml.gz: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.xml.gz
#uniprot_trembl_plants.xml.gz: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_plants.xml.gz
python ~/script/anno_dir/uniprot_xml2info.py -i uniprot_sprot_plants.xml.gz -o uniprot_sprot_plants.fa -a uniprot_sprot_plants.tsv
diamond makedb --in uniprot_sprot_plants.fa --db uniprot_sprot_plants
python ~/script/anno_dir/uniprot_xml2info.py -i uniprot_trembl_plants.xml.gz -o uniprot_trembl_plants.fa -a uniprot_trembl_plants.tsv
diamond makedb --in uniprot_trembl_plants.fa --db uniprot_trembl_plants

## eggNOG (Using Streptophyta here, not all DataBase)
#Streptophyta (35493): http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/35493/
tar -xvf 35493_raw_algs.tar
gunzip -c ./35493/*.gz >> Streptophyta.fa
sed -i 's/-//g' Streptophyta.fa
diamond makedb --in Streptophyta.fa --db Streptophyta
gunzip 35493_annotations.tsv.gz
awk 'BEGIN{OFS="\t"} {print $2, "[" $3 "]"" " $4}' 35493_annotations.tsv > og2anno.tsv
gunzip -c 35493_members.tsv.gz | cut -f2,5 | awk -F'\t' '{n=split($2,a,",");for(i=1;i<=n;i++)print $1"\t"a[i]}' > id2og.tsv
awk -F'\t' 'NR==FNR{a[$2]="["$3"] "$4;next} $1 in a{print $2"\t"a[$1]}' 35493_annotations.tsv id2og.tsv > id2anno.tsv
tar -zxvf 35493_hmms.tar.gz
for file in 35493/*.hmm; do filename=$(basename ${file} .hmm); sed -i "3i ACC   $filename\nDESC  $filename" "$file"; done
ls 35493/*.hmm |  serialize_hmm_models -b Streptophyta.pkl.gz
rm -rf 35493