#!/bin/bash

pep=$1
PWD=$(pwd)
seqkit seq -n $pep | awk '{print $1}' | sed '1i\gene' > allpep.id
sed -i 's/*//g' $pep
## NR
diamond blastp -d ~/00.databases/ncbi/NR_Magnoliopsida.dmnd -q $pep --sensitive --evalue 1e-10 --out nr.blast.out
python -m jcvi.formats.blast best -n 1 nr.blast.out
awk -F'\t' 'NR==FNR{a[$1]=$2;next}{if($2 in a) print $1"\t"a[$2]}' ~/00.databases/ncbi/NR_Magnoliopsida_id2anno.tsv nr.blast.out.best | sed '1i\gene\tNR' > NR.anno
rm nr.blast.out.best
## uniprot_sport_plant
diamond blastp -d ~/00.databases/uniprot/uniprot_sprot_plants.dmnd -q $pep --ultra-sensitive --evalue 1e-10 --out Swissprot.blast.out
python -m jcvi.formats.blast best -n 2 Swissprot.blast.out
awk '{a[$1]=a[$1] (a[$1]?",":"") $2} END{for(i in a) print i"\t"a[i]}' Swissprot.blast.out.best > id2usp.tsv
python ~/script/anno_dir/uniprot2anno.py -i id2usp.tsv -a ~/00.databases/uniprot/uniprot_sprot_plants.tsv -o Swissprot.anno
rm Swissprot.blast.out.best id2usp.tsv
## uniprot_trembl_plants
diamond blastp -d ~/00.databases/uniprot/uniprot_trembl_plants.dmnd -q $pep --sensitive --evalue 1e-10 --out TrEMBL.blast.out
python -m jcvi.formats.blast best -n 2 TrEMBL.blast.out
awk '{a[$1]=a[$1] (a[$1]?",":"") $2} END{for(i in a) print i"\t"a[i]}' TrEMBL.blast.out.best > id2upt.tsv
python ~/script/anno_dir/uniprot2anno.py -i id2upt.tsv -a ~/00.databases/uniprot/uniprot_trembl_plants.tsv -o TrEMBL.anno
rm TrEMBL.blast.out.best id2upt.tsv
## KEGG & Pathway
exec_annotation --cpu 128 -f detail-tsv -p ~/00.databases/kofam/profiles/ -k ~/00.databases/kofam/ko_list -E 1e-5 -o kegg.txt --tmp-dir ./kegg_tem $pep 
# "pykofamsearch" can be used (https://github.com/jolespin/pykofamsearch)
python ~/script/anno_dir/post-kofamscan.py -i kegg.txt -k ~/script/anno_dir/Richer_anno/ko2pathway.list -d ~/script/anno_dir/Richer_anno/pathway_descript.txt -o KEGG.anno -b gene.pathway.bg
## eggNOG
pyhmmsearch -i $pep -b ~/00.databases/emapper/Streptophyta.pkl.gz -o eggnog.hmm.out -p 128 -e 1e-10
reformat_pyhmmsearch -b -i eggnog.hmm.out -o eggnog.hmm.out.best
awk -F'\t' 'NR==FNR{a[$1]=$2;next} FNR==1{next} {if($2 in a) print $1"\t"a[$2]}' ~/00.databases/emapper/og2anno.tsv eggnog.hmm.out.best | sed '1i\gene\teggNOG' > eggNOG.anno
#diamond blastp -d ~/00.databases/emapper/Streptophyta.dmnd -q $pep --sensitive --evalue 1e-10 --out eggnog.blast.out
#python -m jcvi.formats.blast best -n 1 eggnog.blast.out
#awk -F'\t' 'NR==FNR{a[$1]=$2;next}{if($2 in a) print $1"\t"a[$2]}' ~/00.databases/emapper/id2anno.tsv eggnog.blast.out.best | sed '1i\gene\teggNOG' > eggNOG.anno
#rm eggnog.blast.out.best
## InterproScan
nohup ~/software/interproscan-5.70-102.0/interproscan.sh -cpu 128 -f tsv -dp -goterms -pa -iprlookup -i $pep -o interproscan &
# GO
python ~/script/anno_dir/interpro2GO.py -i interproscan -a ~/script/anno_dir/Richer_anno/GOid2altid.txt -d ~/script/anno_dir/Richer_anno/GO_descript.txt -o GO.anno -g gene.GOrich.bg
# interpro_accession
cut -f1,12,13 interproscan | grep "IPR" | awk -F'\t' '{print $1"\t"$2":"$3}' | sort -k1.9n | uniq | awk -F"\t" '{a[$1]=a[$1]$2"; "} END {for (i in a) {print i"\t"a[i]}}' | sed 's/; $//g' | sort -k1.9n | sed '1i\gene\tInterpro' > INTERPRO.anno
# pfam
awk -F'\t' '$4 == "Pfam" {print $1"\t"$5":"$6}' interproscan | sort -k1.9n | uniq | awk -F"\t" '{a[$1] = (a[$1] ? a[$1] "; " : "") $2} END {for (i in a) print i "\t" a[i]}' | sort -k1.9n | sed '1i\gene\tPFAM' > PFAM.anno
# PANTHER
awk -F'\t' '$4 == "PANTHER" {print $1"\t"$5":"$6}' interproscan | sort -k1.9n | uniq | awk -F"\t" '{a[$1] = (a[$1] ? a[$1] ";" : "") $2} END {for (i in a) {sub(/; $/, "", a[i]); print i "\t" a[i]}}' | sort -k1.9n | sed '1i\gene\tPATHER' > PATHER.anno
## iTAK
cd ~/software/iTAK/
pixi run itak -- $PWD/$pep
cd $PWD
awk -F'\t' '{print $1"\t"$2" ["$3"]"}' ${pep}_output/tf_classification.txt | sed '1i\gene\tiTAK_TF' > iTAK_TF.anno
cut -f1,2 ${pep}_output/shiu_classification.txt | sed '1i\gene\tiTAK_PK' > iTAK_PK.anno
## INFO
total=$(awk 'NR>1 && $1!=""{print $1}' allpep.id | sort -u | wc -l)
echo -e "Database\tAnnotated\tRate" > annotation.summary
for ann in *.anno; do
    n=$(awk 'NR>1 && $1!=""{print $1}' "$ann" | sort -u | wc -l)
    rate=$(awk -v n="$n" -v t="$total" 'BEGIN{printf "%.2f%%",n/t*100}')
    echo -e "${ann%.anno}\t$n\t$rate" >> annotation.summary
done
all=$(for ann in *.anno; do awk 'NR>1{print $1}' "$ann"; done|sort -u|wc -l); 
awk -v n=$all -v t=$total 'BEGIN{printf "ALL\t%d\t%.2f%%\n",n,n/t*100}' >> annotation.summary
## combine
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' NR.anno allpep.id > all.NR.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' GO.anno allpep.id | cut -f2 > all.GO.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; c[$1]=$3; next} {if ($1 in b) {$2=b[$1]; $3=c[$1]} else {$2="-"; $3="-"}; print}' OFS='\t' KEGG.anno allpep.id | cut -f2,3 > all.KEGG.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' eggNOG.anno allpep.id | cut -f2 > all.eggNOG.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' INTERPRO.anno allpep.id | cut -f2 > all.INTERPRO.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' PATHER.anno allpep.id | cut -f2 > all.PATHER.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' PFAM.anno allpep.id | cut -f2 > all.PFAM.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' iTAK_TF.anno allpep.id | cut -f2 > all.iTAK_TF.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' iTAK_PK.anno allpep.id | cut -f2 > all.iTAK_PK.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' Swissprot.anno allpep.id | cut -f2 > all.Swissprot.anno
awk -F'\t' 'NR==FNR {b[$1]=$2; next} {if ($1 in b) $2=b[$1]; else $2="-"; print}' OFS='\t' TrEMBL.anno allpep.id | cut -f2 > all.TrEMBL.anno
paste -d "\t" all.NR.anno all.GO.anno all.KEGG.anno all.PFAM.anno all.iTAK_TF.anno all.iTAK_PK.anno all.eggNOG.anno all.PATHER.anno all.INTERPRO.anno all.Swissprot.anno all.TrEMBL.anno > annotation.tsv
rm all.*
