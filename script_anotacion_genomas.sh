#!/bin/bash

#Determine ORFs using prodigal to perform functional domain annotation
for sequence in genomes/*.fa*; do (sequence_name="${sequence//*"/"/}" && sequence_name="${sequence_name//.*/}" && prodigal -f gff -a "${sequence_name}"_aaORFs.fasta -i $sequence -o "${sequence_name}"_ORFs.gff -d "${sequence_name}"_ORFs.gff); done
#Store prodigal results in different folders
mkdir prodigal_results
mv *ORFs* prodigal_results/

#Annotate genomes - general functional domain annotation
parallel -j 20 'hmmsearch --cpu 1 -E 1e-10 --tblout {/.}_pfam_annotation_tbl -o {/.}_pfam_annotation /database/dir/Pfam-A.hmm {}' ::: prodigal_results/*_aaORFs.fasta
#Store hmmer outputs in a different folder
mkdir hmmer_results
mkdir hmmer_results/annotation
mkdir hmmer_results/tbl
mv *_pfam_annotation hmmer_results/annotation/
mv *_pfam_annotation_tbl hmmer_results/tbl/

#Annotate genomes - carbohydrate active enzymes
parallel -j 20 'run_dbcan {} --db_dir /database/dir/cazy_db prok --out_dir {/.}' ::: genomes/*.fa*
#Store cazy outputs in a different folder
mkdir cazy_results
mv *_GCA_* cazy_results/
#Create folder to store functional domains results in csv
mkdir cazy_functional_domains
for folder in cazy_results/*/ ; do (cd "$folder" && code=$(basename "$PWD") && mv hmmer.out "${code}"_hmmer.csv); done
cp cazy_results/*/*_hmmer.csv cazy_functional_domains
#Check all files were correctly copied
ls cazy_functional_domains/*_hmmer.csv | wc -l




#Finish
