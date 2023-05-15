#!/bin/bash

##################################################################

### CONDA ENVIRONMENT ACTIVATION

eval "$(conda shell.bash hook)" # Allows for conda environment use from script

# Check if "chemoreceptor_pipeline" conda environment is present
conda env list > list_env_conda.txt
{ if ! grep -q "chemoreceptor_pipeline" list_env_conda.txt ; then echo "### ERROR : Please install the conda environment (see the README.md), EXITING SCRIPT" ; rm -f list_env_conda.txt ; exit 0 ; fi ; }
rm -f list_env_conda.txt

conda activate chemoreceptor_pipeline

##################################################################

### ARGUMENTS AND VARIABLES

# Run date
run_date=$(date '+%Y%m%d')

# Arguments (with spaces substitution in the database, alignement and genome files -> no spaces in sequence names)
genome=$1 
target_database_filtered=$2 
target_database_full=$3 
extended_database=$4 
align_verif_file=$5 
scripts_folder_location=$6 ; scripts_location=$( echo "$scripts_folder_location" | sed 's/\/$//' )
maximum_intron_length=$7 
number_of_threads=$8
tm_prediction=$9
family=${10} # target gene family, CASE SENSITIVE, used in greps to distinguish family gene from outgroups
exonerate_exhaustive=${11}

# In-script variables
blast_Evalue="1e-4"
getorf_minsize="1000"

# Infos on genome
echo """

################## ################## ################## ################## 
################## ################## ################## ################## 
$(date '+%d/%m/%Y %H:%M:%S')
Genome : "$genome"
START
################## ################## ################## ################## 
################## ################## ################## ################## 

################## Arguments given: ##################
genome file =              $genome
target database filtered = $target_database_filtered
target database full =     $target_database_full
extended database =        $extended_database
alignment outgroup =       $align_verif_file
scripts folder location =  $scripts_folder_location
maximum intron length =    $maximum_intron_length
number of threads =        $number_of_threads
7tm domain filter =        $tm_prediction
target gene family =       $family
"""

##################################################################

### Test for correct argument values before starting the pipeline, else exits (in {} to allow script exit from within if condition)
{ if ! test -e "$genome" ; then echo "### ERROR : genome file doesn't exist, EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! test -e "$target_database_filtered" ; then echo "### ERROR : family database (filtered) file doesn't exist, EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! test -e "$target_database_full" ; then echo "### ERROR : family database (full) file doesn't exist, EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! test -e "$extended_database" ; then echo "### ERROR : extended database file doesn't exist, EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! test -e "$align_verif_file" ; then echo "### ERROR : alignment outgroup file doesn't exist, EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! test -d "$scripts_folder_location" ; then echo "### ERROR : wrong scripts folder name or localization, EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! [[ $maximum_intron_length =~ ^[0-9]+$ ]] ; then echo "### ERROR : maximum intron length must be an integer, EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! [[ $number_of_threads =~ ^[0-9]+$ ]] ; then echo "### ERROR : number of threads must be an integer, EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! ( [[ $tm_prediction == "True" ]] || [[ $tm_prediction == "False" ]] ) ; then echo "### ERROR : DeepTMHMM filter argument must be either 'True' or 'False', EXITING SCRIPT" ; exit 0 ; fi ; }
{ if test -z "$family" ; then echo "### ERROR : please input the targetted family identifier (i.e. 'OR-Receptor'), EXITING SCRIPT" ; exit 0 ; fi ; }
{ if ! grep -q "$family" "$target_database_full"  ; then echo "### ERROR : the targetted family identifier is not found in the target family database, EXITING SCRIPT" ; exit 0 ; fi ; }
### Test for execution rights, with Tree_parser.R as a test sample
{ if ! test -x "$scripts_folder_location"/Tree_parser.R ; then echo "### ERROR : Please run 'chmod u+x -R *' to authorize scripts execution in the pipeline directory and sub-directories" ; exit 0 ; fi ; }

##################################################################

### replace spaces, commas, : by underscores (except in alignment file)
files_array=( "$genome" "$target_database_filtered" "$target_database_full" "$extended_database" )
for file in ${files_array[@]} ; do
	sed 's/ /_/g' "$file" | sed 's/,/_/g' | sed 's/:/_/g' > temp.file ; mv temp.file "$file"
done

### replace - by underscores in genome only (except in alignment file)
sed 's/-/_/g' "$genome" > temp.file ; mv temp.file "$genome"

rm -f temp.file # remove temporary file

##################################################################

#makeblastdb so we can blast genes against the genome ( + samtools faidx index )
if test -f "$genome".ndb ; then echo "Genome blast database index already exists" ; else makeblastdb -in "$genome" -dbtype nucl ; fi 
if test -f "$genome".fai ; then echo "Genome fai file already exists" ; else samtools faidx "$genome" ; fi 

#makeblastdb of target family (full/filtered) and target+outgroup db if needed
if test -f "$target_database_filtered".pdb ; then echo "Target (filtered) blast database index already exists" ; else makeblastdb -in "$target_database_filtered" -dbtype prot ; fi 
if test -f "$target_database_full".pdb ; then echo "Target (full) blast database index already exists" ; else makeblastdb -in "$target_database_full" -dbtype prot ; fi 
if test -f "$extended_database".pdb ; then echo "Target + outgroup blast database index already exists" ; else makeblastdb -in "$extended_database" -dbtype prot ; fi 


#Perform tblastn using known target genes against the genome with a set evalue
tblastn -query "$target_database_filtered" -db "$genome" -evalue $blast_Evalue -outfmt 6 -out target_vs_Genome.blastn -num_threads $number_of_threads


#Lets launch a Rscript that will merge all blast hits. Results file : Blast_nonoverlapping.tsv 
Rscript "$scripts_location"/Rscript_merge_blast_hits.R

xargs samtools faidx "$genome" < Blast_nonoverlapping.tsv > Blast_nonoverlapping.fasta


#Retain only best hits that best match to a target gene
blastx -query Blast_nonoverlapping.fasta -db "$extended_database" -max_target_seqs 1 -outfmt '6 qseqid sseqid' -out blastx_blast_regions.tsv -num_threads $number_of_threads &> /dev/null
grep "$family" blastx_blast_regions.tsv | cut -f1 | sort -u > target_best_hits.txt
rm -f target_Regions.tsv
for i in $( cat target_best_hits.txt ) ; do grep "$i" Blast_nonoverlapping.tsv >> target_Regions.tsv ; done


#Extend all best hits by upstream and downstream . Result file : Potential_target_regions.tsv
Rscript "$scripts_location"/Rscript_merge_filter_extend_blast_hit.R $maximum_intron_length


#Split the target database and launch exonerate with these sequences against potential target regions
mkdir -p Exonerate_split_db
"$scripts_location"/exonerate-2.2.0-x86_64/bin/fastasplit -f "$target_database_filtered" -c $number_of_threads --output Exonerate_split_db


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Search for target family genes in a loop  ###################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
echo """

################## ################## ##################
################## Search for $family genes in a loop 
################## ################## ##################

"""


#re-initialize files

rm -f Potential_multiple_exon_CDS.fa 
rm -f Pseudogenes_multiple_exon.fa 
rm -f No_target_genes_coordinates.txt
rm -f Frameshift_less_Pseudogenes.fa

#Start the loop to search for target genes
loop_nb=0
current_nb_sequences=1
previous_iteration_nb_sequences=$( if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi )
number_regions_blast=$( grep "[0-9]" Potential_target_regions.tsv | wc -l )


#while [ "$current_nb_sequences" -gt "$previous_iteration_nb_sequences" ] && [ "$number_regions_blast" -gt "0" ] ; do
while [ "$current_nb_sequences" -gt "$previous_iteration_nb_sequences" ] && [ "$number_regions_blast" -gt "0" ] ; do

	loop_nb=$(( $loop_nb + 1))
	echo """
################## Loop round number $loop_nb ##################
	"""

	previous_iteration_nb_sequences=$( if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ;  fi )
	
	#Extract identified regions in a fasta file
	xargs samtools faidx "$genome" < Potential_target_regions.tsv > Potential_target_regions.fa
	
	echo "Exonerate running -- May take a while"

	rm -f Exonerate_results.txt
	{ ls Exonerate_split_db/* | parallel -j$number_of_threads "$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate $exonerate_exhaustive --verbose 0 --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' $( echo {} | sed 's/Exonerate_split_db\///g' ) Potential_target_regions.fa >> Exonerate_results.txt
	} &> /dev/null
	sleep 2 # necessary for the calculus servers not to continue without all exonerate chunks done

	wait # necessary not to continue without waiting for every exonerate chunk to be done

	echo "Exonerate research done"

	#extract vulgar lines
	grep "vulgar" Exonerate_results.txt > vulgar_lines.txt 
	
	
	#extract interesting columns of vulgar lines
	#query, query_start, query_end, scaffold, scaffold_start, scaffold_end, strand, exonerate_score
	sed 's/vulgar: //g' vulgar_lines.txt | cut -f1,2,3,5,6,7,8,9 -d " " > vulgar_lines_parsed.txt
	
	
	#count the number of introns using vulgar lines
	IFS=$'\n'
	awk -F'|' 'BEGIN{print "count", "lineNum"}{print gsub(/ I /,"") "\t" NR}' vulgar_lines.txt > number_introns_per_line.txt
	grep -v "count" number_introns_per_line.txt | cut -f1  > intron_numbers.txt 
	
	#add the intron number to vulgar lines 
	paste -d " " vulgar_lines_parsed.txt intron_numbers.txt > vulgar_lines_intron_numbers.txt
	
	
	##Add informations about the best blastp results of each exonerate predicted genes

	sed -n '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/p' Exonerate_results.txt | sed 's/C4 Alignment:.*//g' | sed 's/Hostname:.*//g' | sed 's/Command line:.*//g' | sed 's/^--$//g' | sed 's/-- completed exonerate analysis.*//g' | sed 's/# --- END OF GFF DUMP ---//g' | sed 's/^#$/>seq_to_rename/g' > List_exonerate_cds.fasta #extract all predicted genes sequences
	transeq List_exonerate_cds.fasta List_exonerate_cds.prot &> /dev/null #translate sequences 
	sed 's/ /_/g' vulgar_lines_intron_numbers.txt > sequences_names.txt #extract exonerate vulgar line to rename sequences
	awk '/^>/ { printf("%s_%s\n",$0,i++);next;} { print $0;}' List_exonerate_cds.prot > List_exonerate_cds_renamed.prot #first round of rename
	grep ">" List_exonerate_cds_renamed.prot | sed 's/>//g' > old_names.txt #extract names
	paste -d "\t" old_names.txt sequences_names.txt > renaming_file #create a file for Rename_fasta.pl
	perl "$scripts_location"/Rename_fasta.pl renaming_file List_exonerate_cds_renamed.prot > List_exonerate_cds.prot #completely rename sequences with the exonerate vulgar line
	
	#Perform the blastp
	blastp -query List_exonerate_cds.prot -db "$target_database_full" -outfmt '6 qseqid sseqid evalue' -out all_blastp.txt -max_target_seqs 1 -num_threads $number_of_threads &> /dev/null
	
	#Extract the information
	rm -f all_blastp_parsed.txt
	for i in $( cat sequences_names.txt ) ; do if grep -q "$i" all_blastp.txt ; then grep -m1 "$i" all_blastp.txt | cut -f2,3 >> all_blastp_parsed.txt ; else echo "NoQuery	99999" >> all_blastp_parsed.txt ; fi ; done 
	sed -i 's/	/ /g' all_blastp_parsed.txt
	paste -d " " vulgar_lines_intron_numbers.txt all_blastp_parsed.txt > vulgar_lines_intron_numbers_blastrslt.txt
	
	
	
	#Parse exonerate results. Find the best exonerate results that are most likely complete genes or pseudogenes, and not overlapping
	#This R script will also remove genes that are most likely the merge of two real genes (removed if there are more tblastn non-overlapping results than the number of exons predicted by exonerate)
	Rscript "$scripts_location"/Parse_exonerate_results.R 	#result file : Parsed_exonerate_gene_regions.tsv
	

	#Parse the R result file ...
	nb_row_parsed_exonerate=$( wc -l < Parsed_exonerate_gene_regions.tsv )
	if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then

		IFS=$'\n'
		
		for line in $( cat Parsed_exonerate_gene_regions.tsv ) ; do
			
			query=$( echo "$line" | cut -f7 )
			scaffold=$( echo "$line" | cut -f1 )
			scaff_start=$( echo "$line" | cut -f2 )
			scaff_end=$( echo "$line" | cut -f3 )
		
			echo "$scaffold:$scaff_start-$scaff_end	$query"
		
		done > Correct_coordinates_for_exonerate.tsv

		#Let's now predict genes more precisely on these regions !
		IFS=$'\n'
		mkdir -p Genes_predictions
		
		for line in $( cat Correct_coordinates_for_exonerate.tsv ) ; do 
			
			scaffold_s_e=$( echo "$line" | cut -f1 )
			best_query=$( echo "$line" | cut -f2 )
			scaffold_s_e_n=$( echo "$line" | cut -f1 | sed 's/:/-/g' )
		
			samtools faidx "$genome" $scaffold_s_e > scaffold.fa
			sed -i 's/:/-/g' scaffold.fa

			samtools faidx "$target_database_full" $best_query > query.prot

			{
			"$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate -E True --verbose 0 --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
			} &> /dev/null

			cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta
		
		
			#extract only the best result if there are two with the same score
			if [ $( grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate ) -ge 2 ] ; then
				sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
				sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
				cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
			fi

		done
		
		
		### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
		
		#Result folder
		mkdir -p Filtered_predictions	
		
		for file in Genes_predictions/*.exonerate ; do

			#extract some infos from file name
			file_name=$( echo "$file" | sed 's/.*\///g' )
			file_name_reduced=$( echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' )
			fasta_file_name=$( echo "$file_name" | sed 's/exonerate/fasta/g' )
			initial_header=$( grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' )
		
			#Test if the predicted gene is a target gene or not
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot &> /dev/null
			blastp -query predicted_cds.prot -db "$extended_database" -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_threads &> /dev/null
			
			#Lets continue only if the best match is an target
			if grep -q "$family" blastp_result ; then
		
				#Define the scaffold  name
				scaffold=$( echo "$file" | sed 's/.*\///g' | sed 's/-.*//g' )
				#Define the strand on which the predicted gene is
				strand=$( grep "	similarity	" $file | cut -f7 )
				#Define the first position of the query on the target sequence
				first_hit_range=$( grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
				#Define the last position of the query on the target sequence
				second_hit_range=$( grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " " )
				
				#Lets extract CDS if the gene is on the negative strand
				if [ $strand == "-" ] ; then 
		
					#If strand is minus, then the first position is:
					target_end=$((first_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((first_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
				
					#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in $( cat target_seq.tsv ) ; do begin_exon=$( echo "$exons" | cut -f2 ) ; end_exon=$( echo "$exons" | cut -f3) ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
		
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $getorf_minsize -find 3 &> /dev/null
					if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=$( grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g' ) ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
				
		
					#Rename the fasta file
					if [ $( grep -c ">" Filtered_predictions/$file_name_reduced.ORF ) -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP &> /dev/null
						
						{
						"$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate -E True --verbose 0 --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						} &> /dev/null

						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else "$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


						extracted_scaffold_start=$( echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2 )
						cds_end_extract=$( grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
						cds_start_extract=$( grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " " ) 
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=$( grep "	exon	" verif_coord.exo | wc -l )
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
		
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon). This is a bit redundant with the Rscript previously used, but sometime I detect merges here and not in the Rscript....
						samtools faidx "$genome" $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl &> /dev/null
						number_blast_hit=$( tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue $blast_Evalue -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l )
						if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
		
		
					#If not ORF found, then determinate the gene state
					elif [ $( grep -c ">" Filtered_predictions/$file_name_reduced.ORF ) -lt 1 ] ; then 
				
						stop_codon_state="FALSE"
						edge_state="FALSE"
						frameshift_state="FALSE"
				
						##Stop codon checking
				
						#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
						transeq predicted_cds.fa predicted_cds.prot &> /dev/null
				
						#Estimate the interval on which we will search stop codons. 
						query_name=$( grep "Query: " $file | sed 's/.*Query: //g' )
						query_total_length=$( grep -m1 "$query_name" "$target_database_full".fai | cut -f2 )
						query_start_position=$( grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
						five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
				
						#Lets see if we find stop codon before the five_percent_position
						stop_codon_nb=$( grep -v ">" predicted_cds.prot | fold -w1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l ) #number of stop codons before the ten percent pos
				
						if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
				
				
						##Frameshift checking
				
						#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
						#We remove border exons if there are less than 60nt in length. Run as iteration.
						grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
						awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
				
						#Check for the presence of frameshift
						frameshift_nb=$( grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}' )
						if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
						if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
				
				
				
						##Edge checking
				
						#Check if the gene is at a conting border
						#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
						gene_start_coord=$( cut -f1 Correct_exons.txt | sort -n | head -1 )
						gene_end_coord=$( cut -f2 Correct_exons.txt | sort -n | tail -1 )
				
						extracted_scaffold_start=$( grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2 )
				
						true_start_coord=$((extracted_scaffold_start + gene_start_coord))
						true_end_coord=$((extracted_scaffold_start + gene_end_coord))
				
						#First check if these coordinates are near the end of scaffolds (<5000 bp)
						if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
						scaffold_length=$( grep -m1 "^$scaffold	" "$genome".fai | cut -f2 ) #extract scaffold length from .fai file
						diff_lengths=$((scaffold_length - true_end_coord))
						if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
						
						#Now check if there are consecutive N near the gene that could indicate conting end
						extanded_start_coord=$((true_start_coord - 200))
						extanded_end_coord=$((true_end_coord + 200))
				
						#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
						consecutive_N_nb=$( samtools faidx "$genome" $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1 )
						if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
						if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
				
				
						##Extract the sequence
				
						rm -f Current_exon_rev.txt
				
						#Extract the corresponding sequence
						for line in $( cat Correct_exons.txt ) ; do

							start_pos=$( echo "$line" | cut -f1 )
							end_pos=$( echo "$line" | cut -f2 )
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
							revseq Current_exon.fa Current_exon_rev.fa
				
							#add the reversed sequence to a text file
							grep -v ">" Current_exon_rev.fa >> Current_exon_rev.txt
				
						done
				
						#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
						exon_nb=$( wc -l Exons_length.txt | sed 's/ .*//g' )
				
						header_name=$( echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g' )
						sed -e "1i>$header_name\\" Current_exon_rev.txt > Filtered_predictions/$file_name_reduced.PSEU
						sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
						cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
						sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
	
	
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx "$genome" $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl &> /dev/null
						number_blast_hit=$( tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue $blast_Evalue -outfmt 6 | awk '{ if ($4 >= 40) { print } }' | wc -l )
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
	
				
					fi
		

				#Lets make the same steps with slight modifications for the + strand
				elif [ $strand == "+" ] ; then 
				
					#If strand is minus, then the first position is:
					target_end=$((second_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((second_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
				
					#Extract the CDS sequence predicted by exonerate and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for exons in $( cat target_seq.tsv ) ; do begin_exon=$( echo "$exons" | cut -f2 ) ; end_exon=$( echo "$exons" | cut -f3 ) ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $getorf_minsize -find 3 -reverse FALSE &> /dev/null
				
				
					#Rename the fasta file 
					if [ $( grep -c ">" Filtered_predictions/$file_name_reduced.ORF ) -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP &> /dev/null

						{
						"$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate -E --verbose 0 --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						} &> /dev/null

						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else "$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


						extracted_scaffold_start=$( echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2 )
						cds_start_extract=$( grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
						cds_end_extract=$( grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " " )
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=$( grep "	exon	" verif_coord.exo | wc -l )
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
					
		
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx "$genome" $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl &> /dev/null
						number_blast_hit=$( tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue $blast_Evalue -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l )
						if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
		
		
					#If not ORF found, then determinate the gene state
				
					elif [ $( grep -c ">" Filtered_predictions/$file_name_reduced.ORF ) -lt 1 ] ; then 
				
						stop_codon_state="FALSE"
						edge_state="FALSE"
						frameshift_state="FALSE"
				
					##Stop codon checking
				
						#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
						transeq predicted_cds.fa predicted_cds.prot &> /dev/null
				
						#Estimate the interval on which we wil search stop codons. 
						query_name=$( grep "Query: " $file | sed 's/.*Query: //g' )
						query_total_length=$( grep -m1 "$query_name" "$target_database_full".fai | cut -f2 )
						query_start_position=$( grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
						five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
				
						#Lets see if we find stop codon before the five_percent_position
						stop_codon_nb=$( grep -v ">" predicted_cds.prot | fold -w1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l ) #number of stop codons before the ten percent pos
				
						if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
				
				
						##Frameshift checking
				
						#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
						#We remove border exons if there are less than 60nt in length. Run as iteration.
				
						grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
						awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
				
						#Check for the presence of frameshift
						frameshift_nb=$( grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}' )
						if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
						if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
				
				
						##Edge checking
				
						#Check if the gene is at a conting border
						#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
						gene_start_coord=$( cut -f1 Correct_exons.txt | sort -n | head -1 )
						gene_end_coord=$( cut -f2 Correct_exons.txt | sort -n | tail -1 )
				
						extracted_scaffold_start=$( grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2 )
				
						true_start_coord=$((extracted_scaffold_start + gene_start_coord))
						true_end_coord=$((extracted_scaffold_start + gene_end_coord))
				
						#First check if these coordinates are near the end of scaffolds (<5000 bp)
						if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
						scaffold_length=$( grep -m1 "^$scaffold	" "$genome".fai | cut -f2 ) #extract scaffold length from .fai file
						diff_lengths=$((scaffold_length - true_end_coord))
						if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
						
						#Now check if there are consecutive N near the gene that could indicate conting end
						extanded_start_coord=$((true_start_coord - 200))
						extanded_end_coord=$((true_end_coord + 200))
				
						#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
						consecutive_N_nb=$( samtools faidx "$genome" $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1 )
						if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
						if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
				
				
						##Extract the sequence
				
						rm -f Current_exon.txt
						#Extract the corresponding sequence
						for line in $( cat Correct_exons.txt ) ; do
						
							start_pos=$( echo "$line" | cut -f1 )
							end_pos=$( echo "$line" | cut -f2 )
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
				
							#add the reversed sequence to a text file
							grep -v ">" Current_exon.fa >> Current_exon.txt
				
						done
				
						#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
						exon_nb=$( wc -l Exons_length.txt | sed 's/ .*//g' )
				
						header_name=$( echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g' )
						sed -e "1i>$header_name\\" Current_exon.txt > Filtered_predictions/$file_name_reduced.PSEU
						sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
						cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
						sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
	
						#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
						samtools faidx "$genome" $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
						makeblastdb -in Verification_scaffold.fa -dbtype nucl &> /dev/null
						number_blast_hit=$( tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue $blast_Evalue -outfmt 6 | awk '{ if ($4 >= 40) { print } }' |  wc -l )
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
						if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi 
	
					fi

				fi
		
			else echo "$initial_header" >> No_target_genes_coordinates.txt
		
			fi
		
		done
	fi
	
	
	#Now that we have filtered all our results, we can concatenate the results
	for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds

	cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
	awk '/^>/{f=!d[$1];d[$1]=1}f' Potential_multiple_exon_CDS.fa > non_redundant.fa ; mv non_redundant.fa Potential_multiple_exon_CDS.fa #removes redundant sequences from Potential_multiple_exon_CDS.fa 

	cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...

	cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa

 
	#Extract coordinates of found genes
	grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
	grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
	if test -f No_target_genes_coordinates.txt ; then sed 's/-/	/g' No_target_genes_coordinates.txt >> Coordinates_already_examined.tsv ; fi
	if [ $( wc -l < Coordinates_already_examined.tsv ) -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_already_examined.tsv ; fi
	
	current_nb_sequences=$( if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ; fi )
	
	#re-process blast result to find potential target regions exclusing already found genes
	Rscript "$scripts_location"/Rscript_merge_filter_extend_blast_hit_Second.R $maximum_intron_length

	rm -rf Filtered_predictions/
	rm -rf Genes_predictions/
	rm Parsed_exonerate_gene_regions.tsv
	
	if test -f "Potential_target_regions.tsv" ; then number_regions_blast=$( grep "[0-9]" Potential_target_regions.tsv | wc -l ) ; else number_regions_blast=0 ; fi
	
	echo """
	Number of potential $family regions after loop $loop_nb : $current_nb_sequences
	Number of potential $family pseudogenes after loop $loop_nb : $( if test -f "Pseudogenes_multiple_exon.fa" ; then grep -c ">" Pseudogenes_multiple_exon.fa ; else echo "0" ; fi ) """

done

cp Coordinates_already_examined.tsv Coordinates_already_examined_after_while_loop.tsv

rm -rf Exonerate_split_db
rm -rf Genes_predictions/
rm -f Parsed_exonerate_gene_regions.tsv


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Search for remaining family genes with length a bit below  ##################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
echo """

################## ################## ################## ##################
################## Search for remaining $family genes with length a bit below
################## ################## ################## ##################

"""

Rscript "$scripts_location"/Parse_exonerate_results_second.R

# Exonerate_parse function. if argument "pseudogenes" given, doesn't execute the code only necessary for the gene part
Exonerate_parse () {
	nb_row_parsed_exonerate=$( wc -l < Parsed_exonerate_gene_regions.tsv )

	if [ "$nb_row_parsed_exonerate" -gt "0" ] ; then

		IFS=$'\n'
		
		for line in $( cat Parsed_exonerate_gene_regions.tsv ) ; do
			
			query=$( echo "$line" | cut -f7 )
			scaffold=$( echo "$line" | cut -f1 )
			scaff_start=$( echo "$line" | cut -f2 )
			scaff_end=$( echo "$line" | cut -f3 )
		
			echo "$scaffold:$scaff_start-$scaff_end	$query"
		
		done > Correct_coordinates_for_exonerate.tsv
		
		
		#Let's now predict genes on these regions !
		IFS=$'\n'
		mkdir -p Genes_predictions
		for line in $( cat Correct_coordinates_for_exonerate.tsv ) ; do 
			
			scaffold_s_e=$( echo "$line" | cut -f1 )
			best_query=$( echo "$line" | cut -f2 )
			scaffold_s_e_n=$( echo "$line" | cut -f1 | sed 's/:/-/g' )

			samtools faidx "$genome" $scaffold_s_e > scaffold.fa
			sed -i 's/:/-/g' scaffold.fa
		
			samtools faidx "$target_database_full" $best_query > query.prot

			{
			"$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate -E True --verbose 0 --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo "%tcs" --bestn 1 query.prot scaffold.fa > Genes_predictions/$scaffold_s_e_n.exonerate
			} &> /dev/null

			cp scaffold.fa Genes_predictions/$scaffold_s_e_n.fasta

			#extract only the best result if there are two with the same score
			if [ $( grep -c "Query: " Genes_predictions/$scaffold_s_e_n.exonerate ) -ge 2 ] ; then
				sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_infos
				sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Genes_predictions/$scaffold_s_e_n.exonerate > first_result_sequence
				cat first_result_infos first_result_sequence > Genes_predictions/$scaffold_s_e_n.exonerate
			fi
		
		done
		
		
		### Now we will extract coding sequences from exonerate files. We will define if predicted genes are functionnal or pseudogene ###
		
		
		#Result folder
		mkdir -p Filtered_predictions	
		
		for file in Genes_predictions/*.exonerate ; do
		
			#extract some infos from file name
			file_name=$( echo "$file" | sed 's/.*\///g' )
			file_name_reduced=$( echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' )
			fasta_file_name=$( echo "$file_name" | sed 's/exonerate/fasta/g' )
			initial_header=$( grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' )
		
			#Test if the predicted gene is a target gene or not
			awk '/# --- END OF GFF DUMP ---/ {p=1}; p; /C4 Alignment:/ {p=0}' $file | grep -v "^-- completed" | grep -v "C4 Align" | grep -v "END OF GFF" | sed "s/#/>predicted_cds/g"  > predicted_cds.fa
			transeq predicted_cds.fa predicted_cds.prot &> /dev/null
			blastp -query predicted_cds.prot -db "$extended_database" -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_threads &> /dev/null
			
		
			#Lets continue only if the best match is a target gene
			#if grep -q -i "olfactory\|odorant" blastp_result ; then 
			if grep -q "$family" blastp_result ; then
		
				#Define the scaffold  name
				scaffold=$( echo "$file" | sed 's/.*\///g' | sed 's/-.*//g' )
				#Define the strand on which the predicted gene is
				strand=$( grep "	similarity	" $file | cut -f7 )
				#Define the first position of the query on the target sequence
				first_hit_range=$( grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
				#Define the last position of the query on the target sequence
				second_hit_range=$( grep -m1 "Target range:" $file | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " " )
				
				#Lets extract CDS if the gene is on the negative strand
				if [ $strand == "-" ] ; then 
				
				#file=Genes_predictions/NC_019879.2-28438421-28440954.exonerate 
		
		
					#If strand is minus, then the first position is:
					target_end=$((first_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((first_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$second_hit_range > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
				
					#Extract the target sequence corresponding to the CDS predicted by exonerate and revseq, and remove fasta header
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for  exons in $( cat target_seq.tsv ) ; do begin_exon=$( echo "$exons" | cut -f2 ) ; end_exon=$( echo "$exons" | cut -f3 ) ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
						
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $getorf_minsize -find 3 &> /dev/null
					if grep -q -i "reverse" Filtered_predictions/$file_name_reduced.ORF ; then sequence_to_grep=$( grep -i "reverse" Filtered_predictions/$file_name_reduced.ORF | sed 's/>//g' | sed 's/ .*//g' )  ; samtools faidx Filtered_predictions/$file_name_reduced.ORF $sequence_to_grep > temporary ; mv temporary Filtered_predictions/$file_name_reduced.ORF ; rm Filtered_predictions/$file_name_reduced.ORF.fai ; else rm Filtered_predictions/$file_name_reduced.ORF ; echo "bad strand" > Filtered_predictions/$file_name_reduced.ORF ; fi
				
					#Rename the fasta file
					if [ $( grep -c ">" Filtered_predictions/$file_name_reduced.ORF ) -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP &> /dev/null

						{
						"$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate -E True --verbose 0 --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						} &> /dev/null

						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else "$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi


						extracted_scaffold_start=$( echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2 )
						cds_end_extract=$( grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
						cds_start_extract=$( grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " " )
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=$( grep "	exon	" verif_coord.exo | wc -l )
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
				
						if [ "$1" != "pseudogenes" ] ; then
							echo "check that there were no merge between two genes with a tblastn"
							#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
							samtools faidx "$genome" $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
							makeblastdb -in Verification_scaffold.fa -dbtype nucl &> /dev/null
							number_blast_hit=$( tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue $blast_Evalue -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l )
							if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
						fi

		
					#If not ORF found, then determinate the gene state
					elif [ $( grep -c ">" Filtered_predictions/$file_name_reduced.ORF ) -lt 1 ] ; then 
				
						stop_codon_state="FALSE"
						edge_state="FALSE"
						frameshift_state="FALSE"
				
						##Stop codon checking
				
						#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
						transeq predicted_cds.fa predicted_cds.prot &> /dev/null
				
						#Estimate the interval on which we wil search stop codons. 
						query_name=$( grep "Query: " $file | sed 's/.*Query: //g' )
						query_total_length=$( grep -m1 "$query_name" "$target_database_full".fai | cut -f2 )
						query_start_position=$( grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
						five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
				
						#Lets see if we find stop codon before the five_percent_position
						stop_codon_nb=$( grep -v ">" predicted_cds.prot | fold -w1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l ) #number of stop codons before the five percent pos
				
						if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
				
				
						##Frameshift checking
				
						#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
						#We remove border exons if there are less than 60nt in length. Run as iteration.
						grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
						awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
				
						#Check for the presence of frameshift
						frameshift_nb=$( grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}' )
						if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
						if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
				

						##Edge checking
				
						#Check if the gene is at a conting border
						#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
						gene_start_coord=$( cut -f1 Correct_exons.txt | sort -n | head -1 )
						gene_end_coord=$( cut -f2 Correct_exons.txt | sort -n | tail -1 )
				
						extracted_scaffold_start=$( grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2 )
				
						true_start_coord=$((extracted_scaffold_start + gene_start_coord))
						true_end_coord=$((extracted_scaffold_start + gene_end_coord))
				
						#First check if these coordinates are near the end of scaffolds (<5000 bp)
						if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
						scaffold_length=$( grep -m1 "^$scaffold	" "$genome".fai | cut -f2 ) #extract scaffold length from .fai file
						diff_lengths=$((scaffold_length - true_end_coord))
						if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
						
						#Now check if there are consecutive N near the gene that could indicate conting end
						extanded_start_coord=$((true_start_coord - 200))
						extanded_end_coord=$((true_end_coord + 200))
				
						#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
						consecutive_N_nb=$( samtools faidx "$genome" $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1 )
						if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
						if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
				
				
						##Extract the sequence
				
						rm -f Current_exon_rev.txt
				
						#Extract the corresponding sequence
						for line in $( cat Correct_exons.txt ) ; do
							start_pos=$( echo "$line" | cut -f1 )
							end_pos=$( echo "$line" | cut -f2 )
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
							revseq Current_exon.fa Current_exon_rev.fa

							#add the reversed sequence to a text file
							grep -v ">" Current_exon_rev.fa >> Current_exon_rev.txt
				
						done
				
						#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
						exon_nb=$( wc -l Exons_length.txt | sed 's/ .*//g' )
				
						header_name=$( echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g' )
						sed -e "1i>$header_name\\" Current_exon_rev.txt > Filtered_predictions/$file_name_reduced.PSEU
						sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
						cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
						sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP
		
		
						if [ "$1" != "pseudogenes" ] ; then
							#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
							samtools faidx "$genome" $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
							makeblastdb -in Verification_scaffold.fa -dbtype nucl &> /dev/null
							number_blast_hit=$( tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue $blast_Evalue -outfmt 6 | awk '{ if ($4 >= 40) { print } }' | wc -l )
							if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
							if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi
						fi
				
					fi


				#Lets make the same steps with slight modifications for the + strand
				elif [ $strand == "+" ] ; then 
				
					#If strand is minus, then the first position is:
					target_end=$((second_hit_range + 1))
					#And we will went to extend this by 500bp to be sure to have the potentiel start codon
					target_extanded_end=$((second_hit_range + 500))
				
					#Extract 500bp downstream and extract the whole current scaffold fasta file untill the start
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:1-$first_hit_range > Extend_three_prime.fa
					samtools faidx Genes_predictions/$fasta_file_name $initial_header:$target_end-$target_extanded_end > Extend_five_prime.fa
				
					#remove fasta header of extanded region files
					grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
					grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
				
					#Extract the CDS sequence predicted by exonerate and remove fasta header
		
					grep "	exon	"  $file | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
					for  exons in $( cat target_seq.tsv ) ; do begin_exon=$( echo "$exons" | cut -f2 ) ; end_exon=$( echo "$exons" | cut -f3 ) ; samtools faidx Genes_predictions/$fasta_file_name $initial_header:$begin_exon-$end_exon >> Correct_cds.fa ; done
					grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
				
					#Merge the three regions files, add a fasta header and then search for an ORF with the same parameters we used for single exon genes
					cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
					sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
					getorf -sequence Complete_extanded_sequence.fa -outseq Filtered_predictions/$file_name_reduced.ORF -minsize $getorf_minsize -find 3 -reverse FALSE &> /dev/null
				
				
					#Rename the fasta file 
					if [ $( grep -c ">" Filtered_predictions/$file_name_reduced.ORF ) -ge 1 ] ; then
				
						transeq Filtered_predictions/$file_name_reduced.ORF Filtered_predictions/$file_name_reduced.ORFP &> /dev/null

						{
						"$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate -E True --verbose 0 --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo
						} &> /dev/null

						if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else "$scripts_location"/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Filtered_predictions/$file_name_reduced.ORFP Genes_predictions/$fasta_file_name > verif_coord.exo ; fi

						extracted_scaffold_start=$( echo "$file" | sed 's/.*\///g' | sed 's/.exonerate//g' | sed 's/-/	/g' | cut -f2 )
						cds_start_extract=$( grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
						cds_end_extract=$( grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " " )
						cds_coord_start=$((extracted_scaffold_start + cds_start_extract))
						cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1))
						exon_number=$( grep "	exon	" verif_coord.exo | wc -l )
						sed -i "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end---$exon_number\_exons/g" Filtered_predictions/$file_name_reduced.ORF
					
						if [ "$1" != "pseudogenes" ] ; then

							#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
							samtools faidx "$genome" $scaffold:$cds_coord_start-$cds_coord_end > Verification_scaffold.fa
							makeblastdb -in Verification_scaffold.fa -dbtype nucl &> /dev/null
							number_blast_hit=$( tblastn -query Filtered_predictions/$file_name_reduced.ORFP -db Verification_scaffold.fa -evalue $blast_Evalue -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l )
							if [ "$number_blast_hit" -gt "$exon_number" ] ; then rm Filtered_predictions/$file_name_reduced.ORF ; fi 
						fi
		
		
					#If not ORF found, then determinate the gene state
					elif [ $( grep -c ">" Filtered_predictions/$file_name_reduced.ORF ) -lt 1 ] ; then 
				
						stop_codon_state="FALSE"
						edge_state="FALSE"
						frameshift_state="FALSE"
				
					##Stop codon checking
				
						#lets check for the presence of premature stop codons. We will count the number of predicted stop 5percent before the true end position of the query
						transeq predicted_cds.fa predicted_cds.prot &> /dev/null
				
						#Estimate the interval on which we wil search stop codons. 
						query_name=$( grep "Query: " $file | sed 's/.*Query: //g' )
						query_total_length=$( grep -m1 "$query_name" "$target_database_full".fai | cut -f2 )
						query_start_position=$( grep "Query range: " $file | sed 's/^ *//g' | sed 's/Query range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " " )
						five_percent_position=$((query_total_length * 95 / 100 - query_start_position))
				
						#Lets see if we find stop codon before the five_percent_position
						stop_codon_nb=$( grep -v ">" predicted_cds.prot | fold -w1 | grep -n "\*" | sed 's/:.*//g' | awk -v myvar=$five_percent_position 'BEGIN{FS="\t";OFS="\t"}($1<=myvar){print $1}' | wc -l ) #number of stop codons before the ten percent pos
				
						if [ "$stop_codon_nb" -ge '1' ] ; then stop_codon_state="TRUE" ; fi
				
				
						##Frameshift checking
				
						#To search for frameshifts, start by removing spurious exons at the border. They are most probably true for functionnal genes, and most of the time bad for pseudogenes
						#We remove border exons if there are less than 60nt in length. Run as iteration.
				
						grep "	exon	" $file | cut -f4,5,9 | awk 'BEGIN{FS="\t";OFS="\t"}{{$4=$2-$1} print; }' > Exons_length.txt
						awk 'BEGIN{FS="\t";OFS="\t"}($4>60){print;}' Exons_length.txt > Correct_exons.txt
				
						#Check for the presence of frameshift
						frameshift_nb=$( grep -o "frameshifts [0-9]*" Correct_exons.txt | cut -f2 -d " " | awk '{ sum+=$1} END {print sum}' )
						if [[ $frameshift_nb == "" ]] ; then frameshift_nb=0 ; fi
						if [ "$frameshift_nb" -ge '1' ] ; then frameshift_state="TRUE" ; fi
				
				
						##Edge checking
				
						#Check if the gene is at a conting border
						#These borders are either scaffold end or a repeat of "N", usually more than 50 (100 in zebrafish assembly for example)
						gene_start_coord=$( cut -f1 Correct_exons.txt | sort -n | head -1 )
						gene_end_coord=$( cut -f2 Correct_exons.txt | sort -n | tail -1 )
				
						extracted_scaffold_start=$( grep ">" Genes_predictions/$fasta_file_name | sed 's/>//g' | sed 's/-/	/g' | cut -f2 )
				
						true_start_coord=$((extracted_scaffold_start + gene_start_coord))
						true_end_coord=$((extracted_scaffold_start + gene_end_coord))
				
						#First check if these coordinates are near the end of scaffolds (<5000 bp)
						if [ "$true_start_coord" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
						scaffold_length=$( grep -m1 "^$scaffold	" "$genome".fai | cut -f2 ) #extract scaffold length from .fai file
						diff_lengths=$((scaffold_length - true_end_coord))
						if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
						
						#Now check if there are consecutive N near the gene that could indicate conting end
						extanded_start_coord=$((true_start_coord - 200))
						extanded_end_coord=$((true_end_coord + 200))
				
						#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
						consecutive_N_nb=$( samtools faidx "$genome" $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1 )
						if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
						if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
				
				
						##Extract the sequence
				
						rm -f Current_exon.txt
						#Extract the corresponding sequence
						for line in $( cat Correct_exons.txt ) ; do

							start_pos=$( echo "$line" | cut -f1 )
							end_pos=$( echo "$line" | cut -f2 )
				
							samtools faidx Genes_predictions/$fasta_file_name $initial_header:$start_pos-$end_pos > Current_exon.fa
				
							#add the reversed sequence to a text file
							grep -v ">" Current_exon.fa >> Current_exon.txt
				
						done
				
						#add a header to the text file containing our sequence with the number of exon + frameshift/stopcodon/truncated/edge 
						exon_nb=$( wc -l Exons_length.txt | sed 's/ .*//g' )
				
						header_name=$( echo "$scaffold-$true_start_coord-$true_end_coord---$exon_nb exons-$edge_state-$stop_codon_state-$frameshift_state" | sed 's/ /_/g' )
						sed -e "1i>$header_name\\" Current_exon.txt > Filtered_predictions/$file_name_reduced.PSEU
						sed -i '/^[[:space:]]*$/d' Filtered_predictions/$file_name_reduced.PSEU
						cat predicted_cds.prot > Filtered_predictions/$file_name_reduced.CDSP
						sed -i "s/>.*/>$header_name/g" Filtered_predictions/$file_name_reduced.CDSP

						if [ "$1" != "pseudogenes" ] ; then
							#check that there were no merge between two genes with a tblastn (number of tblastn hits should be inferior or equal to the number of exon)
							samtools faidx "$genome" $scaffold:$true_start_coord-$true_end_coord > Verification_scaffold.fa
							makeblastdb -in Verification_scaffold.fa -dbtype nucl &> /dev/null
							number_blast_hit=$( tblastn -query predicted_cds.prot -db Verification_scaffold.fa -evalue $blast_Evalue -outfmt 6 | awk '{ if ($4 >= 40) { print } }' | wc -l )
							if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.PSEU ; fi 
							if [ "$number_blast_hit" -gt "$exon_nb" ] ; then rm Filtered_predictions/$file_name_reduced.CDSP ; fi
						fi
		
					fi
				fi
		
			else echo "$initial_header" >> No_target_genes_coordinates.txt
		
			fi
		
		done

	fi

	#Now that we have filtered all our results, we can concatenate the results
	for file in Filtered_predictions/*.ORF ; do if grep -q ">" $file ; then i=1 ; else rm $file ; fi ; done #supress files not containg cds
	cat Filtered_predictions/*.ORF >> Potential_multiple_exon_CDS.fa #concatenate results of proper CDS
	cat Filtered_predictions/*.PSEU >> Pseudogenes_multiple_exon.fa #concatenate results of pseudo/edge...
	cat Filtered_predictions/*.CDSP >> Frameshift_less_Pseudogenes.fa

	#Extract coordinates of found genes
	grep ">" Potential_multiple_exon_CDS.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_already_examined.tsv
	grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 >> Coordinates_already_examined.tsv
	if test -f No_target_genes_coordinates.txt ; then sed 's/-/	/g' No_target_genes_coordinates.txt >> Coordinates_already_examined.tsv ; fi
	if [ $( wc -l < Coordinates_already_examined.tsv ) -lt 1 ] ; then echo "Simulated_scaffold	1	10" >> Coordinates_already_examined.tsv ; fi

	current_nb_sequences=$( if test -f "Potential_multiple_exon_CDS.fa" ; then grep -c ">" Potential_multiple_exon_CDS.fa ; else echo "0" ; fi )

	rm -rf Filtered_predictions/
	rm -rf Genes_predictions/
	rm Parsed_exonerate_gene_regions.tsv

}

Exonerate_parse "genes"

cp Coordinates_already_examined.tsv Coordinates_examined_after_smaller_genes_step.tsv

echo """
	Number of potential $family regions after this step : $current_nb_sequences 
	Number of potential $family pseudogenes after this step : $( if test -f "Pseudogenes_multiple_exon.fa" ; then grep -c ">" Pseudogenes_multiple_exon.fa ; else echo "0" ; fi )"""

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Search for target family pseudogenes  #######################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
echo """

################## ################## ##################
################## Search for $family pseudogenes
################## ################## ##################

"""

Rscript "$scripts_location"/Parse_exonerate_results_third.R

Exonerate_parse "pseudogenes"
echo """
	Number of potential $family pseudogene regions after this step : $( if test -f "Pseudogenes_multiple_exon.fa" ; then grep -c ">" Pseudogenes_multiple_exon.fa ; else echo "0" ; fi ) """

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Final round to find target family pseudogenes with length below previous iteration  #########################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
echo """

################## ################## ################## ################## ##################
################## Final round to find $family pseudogenes with length below previous iteration
################## ################## ################## ################## ##################

"""

Rscript "$scripts_location"/Parse_exonerate_results_final.R

Exonerate_parse "pseudogenes"
echo """
	Number of potential $family pseudogene regions after this step : $( if test -f "Pseudogenes_multiple_exon.fa" ; then grep -c ">" Pseudogenes_multiple_exon.fa ; else echo "0" ; fi ) """

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Filter genes and pseudogenes with a blast and a phylogenetic tree  #########################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
echo """

################## ################## ################## ################## ##################
################## Filter $family genes and pseudogenes with a blast and a phylogenetic tree 
################## ################## ################## ################## ##################

"""

#Remove sequences with same headers before 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_multiple_exon.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Pseudogenes_multiple_exon_uniq.fa ; mv Pseudogenes_multiple_exon_uniq.fa Pseudogenes_multiple_exon.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Frameshift_less_Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Frameshift_less_Pseudogenes_uniq.fa ; mv Frameshift_less_Pseudogenes_uniq.fa Frameshift_less_Pseudogenes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Potential_multiple_exon_CDS.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Potential_multiple_exon_CDS_uniq.fa ; mv Potential_multiple_exon_CDS_uniq.fa Potential_multiple_exon_CDS.fa


# We now have all predicted genes in the files Pseudogenes_multiple_exon.fa and Potential_multiple_exon_CDS.fa

#First lets verify our genes with a blast + phylogenetic tree rooted with outgroup genes
grep ">" Pseudogenes_multiple_exon.fa | sed 's/>//g' > Unverified_pseudo_id.txt
fasta_formatter -i Pseudogenes_multiple_exon.fa  > Pseudogenes_multiple_exon_reformat.fa 

#remove sequences with less than 399nt, it does not permit to discriminate below...
seqkit seq -m 133 Frameshift_less_Pseudogenes.fa > Frameshift_less_Pseudogenes_length.fa


#perform a blastx for functional genes and pseudo
blastx -query Potential_multiple_exon_CDS.fa -db "$extended_database" -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out functional_verif_blastx_result -max_target_seqs 1 -num_threads $number_of_threads &> /dev/null
grep "$family" functional_verif_blastx_result | cut -f1 | sort -u > good_blast_complete.txt
xargs samtools faidx Potential_multiple_exon_CDS.fa < good_blast_complete.txt > blast_filtered_Potential_multiple_exon_CDS.fa
mv blast_filtered_Potential_multiple_exon_CDS.fa Potential_multiple_exon_CDS.fa ; rm Potential_multiple_exon_CDS.fa.fai

blastx -query Pseudogenes_multiple_exon.fa -db "$extended_database" -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out verif_blastx_result -max_target_seqs 1 -num_threads $number_of_threads &> /dev/null
grep "$family" verif_blastx_result | cut -f1 | sort -u > good_blast_pseudos.txt
for i in $( cat good_blast_pseudos.txt ) ; do if grep -q "$i" Frameshift_less_Pseudogenes_length.fa ; then echo "$i" >> good_blast_pseudos_length.txt ; fi ; done
xargs samtools faidx Frameshift_less_Pseudogenes_length.fa < good_blast_pseudos_length.txt > good_blast_pseudos_length.fa ; mv good_blast_pseudos_length.fa Frameshift_less_Pseudogenes_length.fa 
rm Frameshift_less_Pseudogenes_length.fa.fai

echo """
################## MAFFT : add predicted pseudogenes and genes to alignment ##################
"""
#Now lets make the tree filter
nb_seq=$( grep -c ">" Frameshift_less_Pseudogenes_length.fa )
if [ "$nb_seq" -gt "0" ] ; then
	mafft --add Frameshift_less_Pseudogenes_length.fa --keeplength "$align_verif_file" > ALL_verification_alignment.aln
else
	cp "$align_verif_file" ./ALL_verification_alignment.aln
fi

#also add functionnal genes
nb_seq=$( grep -c ">" Potential_multiple_exon_CDS.fa )
if [ "$nb_seq" -gt "0" ] ; then
	transeq Potential_multiple_exon_CDS.fa Potential_multiple_exon_CDS.prot &> /dev/null
	sed -i 's/_1$//g' Potential_multiple_exon_CDS.prot
	mafft --add Potential_multiple_exon_CDS.prot --keeplength ALL_verification_alignment.aln > Final_ALL_verification_alignment.aln
else 
	cp ALL_verification_alignment.aln Final_ALL_verification_alignment.aln
fi


#remove sequences with only gaps
python "$scripts_location"/fasta_drop.py Final_ALL_verification_alignment.aln test.fas 1
mv test.fas Final_ALL_verification_alignment.aln

echo """
################## IQTREE : create gene tree with database + predicted genes/pseudogenes ##################
"""

#lets start a tree
iqtree -s Final_ALL_verification_alignment.aln -st AA -nt $number_of_threads -m JTT+F+I+G4 -redo -n 200

### Parse the tree
# known gene (IDs) of the target family
grep ">" $target_database_full | sed 's/>//g' > Known_family_genes_id.txt
# known outgroups (IDs) of the target family
grep ">" $extended_database | grep -v "$family" | sed 's/>//g' > Known_outgroup_genes_id.txt

Rscript "$scripts_location"/Tree_parser.R #This script keep only genes clustering with known target family genes : Current_species_target.txt

grep "FALSE\|TRUE" Current_species_target.txt > target_pseudos.id
grep -v "FALSE\|TRUE" Current_species_target.txt > target_functionnal.id


xargs samtools faidx Potential_multiple_exon_CDS.fa < target_functionnal.id > Functionnal_target_genes.fa
xargs samtools faidx Pseudogenes_multiple_exon_reformat.fa < target_pseudos.id > Pseudogenes_target_genes.fa


## Now lets filter these files : remove ambigous sequences 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_target_genes.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /N/)' | tr "\t" "\n" > clear_Functionnal_target_genes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_target_genes.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /N/)' | tr "\t" "\n" > clear_Pseudogenes_target_genes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Functionnal_target_genes.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /N/)' | tr "\t" "\n" > unclear_Functionnal_target_genes.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_target_genes.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /N/)' | tr "\t" "\n" > unclear_Pseudogenes_target_genes.fa
cat unclear_Functionnal_target_genes.fa unclear_Pseudogenes_target_genes.fa > Ambiguous_target.fa

transeq clear_Functionnal_target_genes.fa clear_Functionnal_target_genes.prot &> /dev/null
sed -i 's/_1$//g' clear_Functionnal_target_genes.prot #translate CDS

if ( [ $tm_prediction == "True" ] || [ $tm_prediction == "TRUE" ] ) ; then
	#Check 7tm domains

	#Check if complete genes have 7tm domain determined with phobius or TMHMM
	echo """
	################## Check 7tm domain : Phobius ##################
	"""
	#First check with phobius
	perl "$scripts_location"/phobius/phobius.pl -long clear_Functionnal_target_genes.prot > Phobius_verification.txt #run phonius in long mode
	grep ">" clear_Functionnal_target_genes.prot | sed 's/>//g' > gene_id.txt #extract cds id

	for gene in $( cat gene_id.txt ) ; do nb_transm=$( sed '/'"$gene"'/,/\/\//!d;/\/\//q' Phobius_verification.txt | grep "TRANSMEM" | wc -l ) ; echo "$gene,$nb_transm" ; done > Gene_NbTm.tsv
	awk 'BEGIN{FS=",";OFS=","}($2>=7){print $1;}' Gene_NbTm.tsv > Phobius_genes_with_7tm.txt
	awk 'BEGIN{FS=",";OFS=","}($2<7){print $1;}' Gene_NbTm.tsv > Phobius_genes_without_7tm.txt


	echo """
	################## Check 7tm domain : DeepTMHMM ##################
	"""
	#Now, with TMHMM --> OLD WAY
	#perl "$scripts_location"/tmhmm-2.0c/bin/tmhmm clear_Functionnal_target_genes.prot > tmhmm_verification.txt
	#grep "Number of predicted TMHs:" tmhmm_verification.txt | sed 's/# //g' | sed 's/ Number of predicted TMHs:  /,/g' | awk -F "," '{ if(($2 >= 7)) { print $1} }' > tmhmm_genes_with_7tm.txt
	#grep "Number of predicted TMHs:" tmhmm_verification.txt | sed 's/# //g' | sed 's/ Number of predicted TMHs:  /,/g' | awk -F "," '{ if(($2 < 7)) { print $1} }' > tmhmm_genes_without_7tm.txt 

	#Now, with DeepTMHMM
	until test -f biolib_results/TMRs.gff3
	do
	biolib run DTU/DeepTMHMM --fasta clear_Functionnal_target_genes.prot 2> /dev/null
	done
	cp biolib_results/TMRs.gff3 ./tmhmm_verification.txt
	grep "Number of predicted TMRs:" tmhmm_verification.txt | sed 's/# //g' | sed 's/ Number of predicted TMRs: /,/g' | awk -F "," '{ if(($2 >= 7)) { print $1} }' > tmhmm_genes_with_7tm.txt
	grep "Number of predicted TMRs:" tmhmm_verification.txt | sed 's/# //g' | sed 's/ Number of predicted TMRs: /,/g' | awk -F "," '{ if(($2 < 7)) { print $1} }' > tmhmm_genes_without_7tm.txt 


	#Combine results of predictions
	cat Phobius_genes_with_7tm.txt tmhmm_genes_with_7tm.txt | sort -u > genes_with_7tm.txt
	sort gene_id.txt > gene_id_sorted.txt ; sort genes_with_7tm.txt > sorted_genes_with_7tm.txt
	comm -3 gene_id_sorted.txt sorted_genes_with_7tm.txt > Genes_without_7tm.txt


	#Print the result in the fasta file
	for gene in $( cat gene_id.txt ) ; do
		pred_phobius="FALSE"
		pred_tmhmm="FALSE"

		if grep -q "$gene" Phobius_genes_with_7tm.txt ; then pred_phobius="TRUE" ; fi
		if grep -q "$gene" tmhmm_genes_with_7tm.txt ; then pred_tmhmm="TRUE" ; fi 

		if [ $pred_phobius == "TRUE" ] && [ $pred_tmhmm == "TRUE" ] ; then 
			new_gene_name="$gene---phobius-deepTMHMM"
		elif [ $pred_phobius == "TRUE" ] && [ $pred_tmhmm == "FALSE" ] ; then 
			new_gene_name="$gene---phobius" 
		elif [ $pred_phobius == "FALSE" ] && [ $pred_tmhmm == "TRUE" ] ; then
			new_gene_name="$gene---deepTMHMM"
		else
			new_gene_name="$gene"
		fi 

		echo "$new_gene_name" 

	done > New_gene_name_with_predictions.txt

	paste -d "\t" gene_id.txt New_gene_name_with_predictions.txt > renaming_file_tm.txt

	perl "$scripts_location"/Rename_fasta.pl renaming_file_tm.txt clear_Functionnal_target_genes.fa > temporary.fasta ; mv temporary.fasta clear_Functionnal_target_genes.fa
fi

# If two genes/pseudogenes are overlapping due to the multiple extensions, then keep only the longest found gene/pseudogene
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Ambiguous_target.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > FINAL_Ambiguous.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Pseudogenes_target_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > clear_Pseudogenes_target_genes_uniq.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Functionnal_target_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > clear_Functionnal_target_genes_uniq.fa

nb_seq=$( grep -c ">" FINAL_Ambiguous.fa )
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" FINAL_Ambiguous.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_ambiguous_final.tsv
fi

nb_seq=$( grep -c ">" clear_Pseudogenes_target_genes_uniq.fa )
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" clear_Pseudogenes_target_genes_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_pseudogenes_final.tsv
fi

nb_seq=$( grep -c ">" clear_Functionnal_target_genes_uniq.fa )
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" clear_Functionnal_target_genes_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_genes_final.tsv
fi

Rscript "$scripts_location"/Remove_redundancy.R


IFS=$'\n'

for line in $( cat best_genes_functionnal.tsv ) ; do scaff=$( echo "$line" | cut -f1 ) ; start=$( echo "$line" | cut -f2 ) ; end=$( echo "$line" | cut -f3 ) ; grep -m1 "$scaff.*$start.*$end" clear_Functionnal_target_genes_uniq.fa | sed 's/>//g' >> functionnal_to_keep.txt ; done
for line in $( cat best_genes_ambigous.tsv ) ; do scaff=$( echo "$line" | cut -f1 ) ; start=$( echo "$line" | cut -f2 ) ; end=$( echo "$line" | cut -f3 ) ; grep -m1 "$scaff.*$start.*$end" FINAL_Ambiguous_target_uniq.fa | sed 's/>//g' >> ambigous_to_keep.txt ; done
for line in $( cat best_genes_pseudogenes.tsv ) ; do scaff=$( echo "$line" | cut -f1 ) ; start=$( echo "$line" | cut -f2 ) ; end=$( echo "$line" | cut -f3 ) ; grep -m1 "$scaff.*$start.*$end" clear_Pseudogenes_target_genes_uniq.fa | sed 's/>//g' >> pseudogenes_to_keep.txt ; done

if test -f "functionnal_to_keep.txt" ; then xargs samtools faidx clear_Functionnal_target_genes_uniq.fa < functionnal_to_keep.txt > clear_Functionnal_target_genes.fa ; else echo "" > clear_Functionnal_target_genes.fa ; fi 
if test -f "ambigous_to_keep.txt" ; then xargs samtools faidx FINAL_Ambiguous_target_uniq.fa < ambigous_to_keep.txt > FINAL_Ambiguous.fa ; else echo "" > FINAL_Ambiguous.fa ; fi 
if test -f "pseudogenes_to_keep.txt" ; then xargs samtools faidx clear_Pseudogenes_target_genes_uniq.fa < pseudogenes_to_keep.txt > clear_Pseudogenes_target_genes.fa ; else echo "" > clear_Pseudogenes_target_genes.fa ; fi 


# copy pseudogenes and genes before (optionnal) filter
cp clear_Pseudogenes_target_genes.fa FINAL_Pseudogenes_no_filter.fa
cp clear_Functionnal_target_genes.fa FINAL_Functionnal_no_filter.fa

#DISCONTINUED UNTIL BETTER PROGRAMS ARE FOUND FOR INSECT OR TM FINDING
#Also classify genes without 7tm as pseudogenes (if the option is TRUE)
#if ( [ $tm_prediction == "True" ] || [ $tm_prediction == "TRUE" ] ) ; then
	#grep ">" clear_Functionnal_target_genes.fa | grep "phobius\|tmhmm" | sed 's/>//g' > 7tm_genes
	#grep ">" clear_Functionnal_target_genes.fa | grep -v "phobius\|tmhmm" | sed 's/>//g' > non_7tm_genes
	#xargs samtools faidx clear_Functionnal_target_genes.fa < 7tm_genes > FINAL_Functionnal_7tm.fa
	#xargs samtools faidx clear_Functionnal_target_genes.fa < non_7tm_genes >> clear_Pseudogenes_target_genes.fa
	#cp clear_Pseudogenes_target_genes.fa FINAL_Pseudogenes_and_no_7tm.fa
#fi


################################################################################################################################################

### prepare result folder name
genome_filename=$( echo "$genome" | sed 's/.*\///' ) 
res_directory=$( echo Results_"$family"_"$genome_filename"_"$run_date" )

### "Stats"
nb_target=$( if test -f "FINAL_Functionnal_no_filter.fa" ; then grep -c ">" FINAL_Functionnal_no_filter.fa ; else echo 0 ; fi )
nb_pseudo=$( if test -f "FINAL_Pseudogenes_no_filter.fa" ; then grep -c ">" FINAL_Pseudogenes_no_filter.fa ; else echo 0 ; fi )
nb_total=$(( $nb_target + $nb_pseudo ))

nb_target_filtered=$( if test -f "FINAL_Functionnal_7tm.fa" ; then grep -c ">" FINAL_Functionnal_7tm.fa ; else echo "NA" ; fi ) 
nb_pseudo_filtered=$( if test -f "FINAL_Pseudogenes_and_no_7tm.fa" ; then grep -c ">" FINAL_Pseudogenes_and_no_7tm.fa ; else echo "NA" ; fi )

nb_ambiguous=$( if test -f "FINAL_Ambiguous.fa" ; then grep -c ">" FINAL_Ambiguous.fa ; else echo 0 ; fi )

nb_truncated=$( if test -f "FINAL_Pseudogenes_no_filter.fa" ; then grep -c "exons-FALSE-FALSE-FALSE" FINAL_Pseudogenes_no_filter.fa ; else echo 0 ; fi )
nb_edges=$( if test -f "FINAL_Pseudogenes_no_filter.fa" ; then grep -c "exons-TRUE" FINAL_Pseudogenes_no_filter.fa ; else echo 0 ; fi )
nb_stop_fs=$(( $nb_pseudo - $nb_edges - $nb_truncated ))

### Add line to the results file

# Create file + header if Results_summary.csv doesn't exist
if ! test -f Results_summary.csv ; then
    echo "Genome,Results_folder,Run_date(YYYYMMDD),Family,Max_intron_length,Total_genes_pseudogenes,Genes,Pseudogenes,Truncated_pseudogenes,Stop_frameshift_pseudogenes,Edge_pseudogenes,Ambiguous" > Results_summary.csv
fi

#Add line for current genome
echo "$genome,$res_directory,$run_date,$family,$maximum_intron_length,$nb_total,$nb_target,$nb_pseudo,$nb_truncated,$nb_stop_fs,$nb_edges,$nb_ambiguous" >> Results_summary.csv

################################################################################################################################################

### END : Infos on genome and results
echo """

################## ################## ################## ################## 
################## ################## ################## ################## 
$(date '+%d/%m/%Y %H:%M:%S')
Genome : "$genome"
END

Number of $family genes :                     $nb_target
Number of $family pseudogenes :               $nb_pseudo """

#DISCONTINUED UNTIL BETTER PROGRAMS ARE FOUND FOR INSECT OR TM FINDING
#if ( [ $tm_prediction == "True" ] || [ $tm_prediction == "TRUE" ] ) ; then
#	echo """
#Number of $family genes with 7tm domain :     $nb_target_filtered ( 7tm filter : ON )
#Number of $family pseudogenes (after filter): $nb_pseudo_filtered ( 7tm filter : ON) """
#fi

echo """
Number of ambiguous genes :                   $nb_ambiguous

Results folder : "$res_directory"/
################## ################## ################## ################## 
################## ################## ################## ################## 

"""

################################################################################################################################################

### Put the result files into results folder for the genome
rm -rf "$res_directory"
mkdir "$res_directory"

if ( [ $tm_prediction == "True" ] || [ $tm_prediction == "TRUE" ] ) ; then
	# rename DeepTMHMM files first
	for filepath in biolib_results/* ; do filename=$( echo "$filepath" | sed 's/.*\///' ) ; mv "$filepath" biolib_results/deepTMHMM_"$filename" ; done
fi

# copy results files in results directory
{ cp biolib_results/* "$res_directory" ; rm -rf biolib_results
cp tmhmm_genes_with_7tm.txt "$res_directory"/deepTMHMM_genes_with_7tm.txt
cp genes_with_7tm.txt "$res_directory"
cp ALL_verification_alignment.aln "$res_directory"
cp Coordinates_* "$res_directory"
cp FINAL_* "$res_directory"
cp Final_ALL_verification_alignment.* "$res_directory"
cp Phobius_genes_with_7tm.txt "$res_directory"
cp Phobius_verification.txt "$res_directory" 
} &> /dev/null

# rename all results files to add genome and target gene family specification
for filepath in  "$res_directory"/* ; do filename=$( echo "$filepath" | sed 's/.*\///' ) ; mv "$filepath" "$res_directory"/"$genome_filename"_"$family"_"$filename" ; done


### Clean all temporary files
#rm -f *.txt *.tsv *.fasta *.fa *.id *.blastn *_genes
#rm -f ALL_* all_* best_* blast* Blast_* clear_* Complete_* Coordinates_* Correct_* Current_* Extend_* Final_* functional_* List_* Potential_* predicted_* Pseudogenes_* renaming_file query.prot Verification_* verif_*


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## END  ########################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

#FINAL (main) RESULTS FILES
#FINAL_Functionnal_no_filter.fa
#FINAL_Functionnal_7tm.fa
#FINAL_Pseudogenes_and_no_7tm.fa
#FINAL_Pseudogenes_no_filter.fa
#FINAL_Ambiguous.fa
