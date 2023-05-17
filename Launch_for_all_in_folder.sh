#!/bin/bash

folder=$( echo "$1" | sed 's/\/$//' )
run_date=$(date '+%Y%m%d')

echo "$folder"

# Test if folder exists
{ if ! test -d "$folder" ; then echo "### ERROR : wrong genome folder name or localization, EXITING SCRIPT" ; exit 0 ; fi ; }

# create file list
file_list=($( ls "$folder" ))

for file in "${file_list[@]}" ; do

    if ! grep -q ">" "$folder/$file" ; then 

        echo "$file is not a nucleotidic fasta, on to next file" 

    else 

        echo "File : $file"
	
       	bash Insect_Chemoreceptor_Finder.sh "$folder/$file" #### ADD THE END OF THE COMMAND HERE; all arguments after "genome_file" (excluded)

        # Delete all genome blast database files
        rm -f "$folder"/*.fna.*
        
    fi 

    echo """
####################################################################################################################################
    """
done

