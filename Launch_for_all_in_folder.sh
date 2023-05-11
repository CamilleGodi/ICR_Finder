#!/bin/bash

folder=$( echo "$1" | sed 's/\/$//' )
run_date=$(date '+%Y%m%d')

echo "$folder"

# Test if folder exists
{ if ! test -d "$folder" ; then echo "### ERROR : wrong genome folder name or localization, EXITING SCRIPT" ; exit 0 ; fi ; }

# create file list
file_list=($( ls "$folder" ))

for file in "${file_list[@]}" ; do

    if ! grep -q ">" "$folder"/"$file" ; then 

        echo "$file is not a nucleotidic fasta, on to next file" 

    else 

        echo "File : $file"
	
       	./ICR_Finder_230510.sh "$folder/$file" DB_Bte_Obi_Nvi/DB_Bte_Obi_Nvi_OR_cdhit70.prot DB_Bte_Obi_Nvi/DB_Bte_Obi_Nvi_OR.prot DB_Bte_Obi_Nvi/DB_Bte_Obi_Nvi_OR_Ame_GR.prot DB_Bte_Obi_Nvi/DB_Bte_Obi_Nvi_OR_Ame_GR.aln Scripts 3000 40 False OR-Receptor -E |& tee "$file"_"$run_date".log

        # Delete all genome blast database files
        rm -f "$folder"/*.fna.*
        
    fi 

    echo """
####################################################################################################################################
    """
done

