# Plan
I/ Introduction  
II/ Intructions  
&nbsp;&nbsp;&nbsp;&nbsp;- A/ Before first use  
&nbsp;&nbsp;&nbsp;&nbsp;- B/ Command details  
&nbsp;&nbsp;&nbsp;&nbsp;- C/ Example  
&nbsp;&nbsp;&nbsp;&nbsp;- D/ Launch all genomes in a folder  
III/ Results  
&nbsp;&nbsp;&nbsp;&nbsp;- A/ Files  
&nbsp;&nbsp;&nbsp;&nbsp;- B/ Interpretation  
VI/ Notes  
&nbsp;&nbsp;&nbsp;&nbsp;- A/ General notes  
&nbsp;&nbsp;&nbsp;&nbsp;- B/ Notes for building custom databases and alignement files  
&nbsp;&nbsp;&nbsp;&nbsp;- C/ Notes on adaptation to different gene or animal families  
V/ Correspondance  
---

---

# I/ Introduction

This pipeline, **Insect chemoreceptors Finder** (ICR_Finder) allows for **automatic and throughout research of chemoreceptors genes/pseudogenes of various families (OR,GR,IR) in insect genomes**. "Database" files are provided for OR research in Hymenoptera genomes.

It can also be adapted to run with custom databases (see part IV/), and supposedly on different animal taxa and gene families (see part V/).

This pipeline is an adaptation and generalisation based on a vertebrates V2R-receptors finder pipeline created by Maxime Policarpo, notably used to generate analysis for *"Coevolution of the olfactory organ and its receptor repertoire in ray-finned fishes"* (BMC Biol., 2022, DOI: 10.1186/s12915-022-01397-x) and other articles. I thank him for the original scripts and for the precious help.

# II/ Intructions

> **IMPORTANT :**   
Do NOT use file names or repertory names containing spaces. (ever !)  
Please avoid to put any .fa file in the working directory, favour the use a data folder. Fasta files are purged from the folder at the end of the pipeline.

## A/ Before first use :
Upon first use, in the pipeline repertory, allow for script execution of all files :
```bash
chmod -R u+x *
```

Then run this script to create the conda environment used by the pipeline :
```bash
conda env create -f conda_chemoreceptor_pipeline.yml
```
> **Note :** if the conda environment is installed on a custom path, you have to provide it to the *conda activate* command in the ICR_Finder.sh script instead of "chemoreceptor_pipeline". If the custom-path-environment isn't in your working directory, you may have to disable the conda environment check in the beginning of the ICR_Finder.sh script.

---

## B/ Command details :

```bash
bash ICR_Finder.sh genome_file path/to/filtered_db path/to/full_db path/to/extended_db path/to/alignment_outgroup path/to/scripts max_intron_length threads_nb tmhmm_filter target_gene_family
```

with:  
| Variable  |  Description |  Format |
|-----------|--------------|---------|
|  **genome_file** | path to the (nucleic) fasta file of your genome | .fa  or .fna file|
| **filtered_db** | path to the chemoreceptor database filtered by cd-hit (representative genes only), amino acid multifasta file | path to .prot, .fa or .faa file
| **full_db** | path to the full chemoreceptor database, amino acid multifasta file | path to .prot, .fa or .faa file
| **extended_db** | path to the full chemoreceptor database with outgroups, amino acid multifasta file | path to .prot or .faa file
| **alignment_outgroup** | path to chemoreceptor-outgroups alignment (see part IV) | .aln file
| **scripts** | path to the Scripts folder | folder path
| **max_intron_length** | maximum intron length tolerated (see part III) | integer |
| **threads_nb** | number of CPU threads used for **Blast**, **Exonerate** and **IQTREE** | integer |
| **7tm_prediction** | if True, adds a transmembrane regions prediction by **Phobius** and **DeepTMHMM**, and associated output files | Boolean : *True* or *False*
| **target_gene_family** | **case sensitive**, term used in the sequence names of genes of the target family (i.e. : OR_receptor) | character string, only alphanumeric characters and underscores/dashes


## C/ Example for Hymenoptera OR prediction

Hymenoptera OR prediction in *genome_file*, with maximum intron size set to 3000, 8 threads used and no 7TM domain prediction :

```bash
bash ICR_Finder.sh ~/Bureau/2_Pipeline_versions/Insect_Chemoreceptors_Finder_Apis_test/Am_no_GPCR/Apis_mellifera.fna DB_HymenopteraORs/DB_HymenopteraORs_ORs_cdhit70.prot DB_HymenopteraORs/DB_HymenopteraORs_ORs.prot DB_HymenopteraORs/DB_HymenopteraORs_ORs_GRs.prot DB_HymenopteraORs/DB_HymenopteraORs_ORs_GRs.aln Scripts 3000 8 False OR-Receptor
```


## D/ Launch all genomes in a folder
To launch the pipeline with the same arguments for a set of genomes contained in a folder (*Folder_name*), substitute the comment of the 24th line in the "Launch_for_all_in_folder.sh" script by the rest of your arguments (see part IIB), then run :

```bash
bash Launch_for_all_in_folder.sh Folder_name
```

# III/ Results
## A/ Files
### General files
All genome results are stored in a folder named *"Results_[target gene family]_[genome file]_[run date]"*, with a prefix reflecting the genome and target gene family.

| File  |  Description |
|-----------|--------------|
| **Results_summary.csv** | table summarizing the results of each run. Unless this file is moved or renamed, new run results are added a new row of the table (In the working directory)
| prefix_**FINAL_Ambiguous.fa** | Ambiguous target family genes (nucleotide fasta) |
| prefix_**FINAL_Functionnal_no_filter.fa** | target family predicted genes (nucleotide fasta) |
| prefix_**FINAL_Pseudogenes_no_filter.fa** | target family predicted pseudogenes (nucleotide fasta) |
| prefix_**Coordinates_already_examined.tsv** | List of all the potential target gene family regions examined by **Exonerate** |
| prefix_**Coordinates_genes_final.tsv** | Predicted target genes coordinates |
| prefix_**Coordinates_pseudogenes_final.tsv** | Predicted target pseudogenes coordinates |
| prefix_**Coordinates_ambiguous_final.tsv** | Predicted target ambiguous genes coordinates |

Files for tree verification : MAFFT alignement results and IQTREE results
| File  |  Description |
|-----------|--------------|
| prefix_**ALL_verification_alignement.aln** | MAFFT alignement file of the predicted proteins with database proteins |
| prefix_**Final_ALL_verification_alignement.aln.\*** | IQTREE diverse files and logs |

### 7tm-filter-related files

If the 7tm (7 transmembrane regions) filter is enabled :

| File  |  Description |
|-----------|--------------|
| prefix_**Phobius_verification.txt** | Phobius transmembrane regions prediction results |
| prefix_**Phobius_genes_with_7tm.txt** | List of the predicted genes with 7tm (or more) according to Phobius |
| prefix_**deepTMHMM_deeptmhmm_results.md** | DeepTMHMM job summary |
| prefix_**deepTMHMM_predicted_topologies.3line** | Detail of the predicted proteins topologies |
| prefix_**deepTMHMM_predicted_TMRs.gff3** | Number and list of predicted proteins transmembrane regions |
| prefix_**deepTMHMM_genes_with_7tm.txt** | List of the predicted genes with 7tm according to DeepTMHMM |
| prefix_**genes_with_7tm.txt** |  List of the predicted genes with 7tm (or more) according to Phobius or DeepTMHMM |

## B/ Interpretation
### Genes
If the 7-transmembrane-regions parameter (**7tm_prediction**) is active, predicted gene names get an additionnal indication if they have 7 transmembrane regions (or more) according to one or both of the protein topology programs used : "phobius", "deepTMHMM".

### Pseudogenes
Each predicted pseudogene has 3 booleans on its name:
- First : **Edge status**. If TRUE : the predicted pseudogene is near to a contig edge (<5000 bp), and there is more than 50 adjacent N nucleotides 200pb before or after this pseudogene.
- Second : **Stop codon status**. If TRUE : the predicted pseudogene has a stop codon outside of the last 5% of its length.
- Third : **Frameshift status**. If TRUE : the predicted pseudogene has 1 or more detected frameshift

Predicted pseudogenes that are FALSE for the three status are categorized as **truncated** in the results csv.

Stop and frameshift predicted pseudogenes in predicted pseudogenes close to an edge (TRUE-\*-\*) are not counted in the *Results_summary.csv* as stop/frameshift as edge genes confuses the **Exonerate** program.

### Ambiguous genes
Predicted genes are categorized as "Ambiguous" when they contain one or more N in their nucleotide sequence.

# IV/ Notes
## A/ General notes
- **max_intron_length** has to be adapted to your dataset to limit overlooking genes with long introns while not merging two close genes together.
- Exonerate is used to identify potential chemoreceptor regions. It may take VERY long!
- DeepTMHMM is currently (as of may 2023) not reliable for 7TM domain prediction.
- If your genome file contains spaces, commas, dashes and/or colons in the sequence names, this pipeline will automatically replace them by underscores.

## B/ Notes for building custom databases and alignement files
- Blast e-value ( in-script variable : "eValueBlast") should be low for large species group , higher for smaller groups, i.e. 1e-4 for Insects and 1e-10 or more for Apis
- For **alignment_outgroup** : clean alignment of the target family database proteins with outgroups proteins. Alignment should be as clean (outgroups well separated from the target family) as possible and representative proteins (after a cd-hit 70/80 filter for example). This file can be obtained through mafft for example (recommanded options : mafft --localpair --maxiterate 2000 ).
- If your genome file contains dashes in the sequence names, this script will automatically replace them by underscores.
- If one of your database files contains spaces and/or colons in the sequence names, this script will automatically replace them by underscores.
- Sequence names should be as unambiguous as possible, with the gene family clearly stated in a distinct manner (ex: OR-Receptor instead of Or/OR), and mention of the species. Use only alphanumeric characters and underscores/dashes.

## C/ Notes on adaptation to different gene or animal families  
> The pipeline was only adapted and tested on Hymenoptera OR chemoreceptors, but with the subsequent modifications, it may work for different gene families or animal taxa.

To change :

### **getorf** ( in-script variable of ICR_Finder.sh )
commands option "*-minsize*" should be adapted to the minimum size of your target protein families, in (exonic) nucleotides size

### ***Parse_exonerate_results* R scripts "*putative protein length*" filters** 
The values need to be tailored to your data. The highest threshold should be slightly lower than the longest proteins known in your target gene family, then decreasing progressively to the minimal values.
The values and repartition should be found by using the pipeline on a model genome with already known repertory for the target gene family, if possible, to see if most genes are found back by the pipeline when using a set of threshold. 

# V/ Correspondance
- Camille GODI, Rouen Normandy University (student, Master Degree Internship), camille.godi[at]univ-rouen.fr or camille.godi[at]gmail.com
- Maxime POLICARPO, Basel University, maxime.policarpo[at]unibas.ch
