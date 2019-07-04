#!/bin/bash          

# Clementine Francois, Faustine Durand
# June 2019
# contact: clementine.francois@univ-lyon1.fr

#set -e


##################################################################################################################################################################################################
######## WHAT THIS SCRIPT DOES :
# for a given species / genome assembly: outputs a list of the contaminant cds and scaffolds (only contamination from non-metazoa). Also outputs a list of potential HGT candidates (to be validated afterwards).
# foreign cds are classified into 5 groups: plant / eubacteria / archaea / fungi / protist
# this script performs a 1st similarity-based (blast) step then a 2nd synteny-based (Gmap) step.
# !! you have to specify the paths to some files: see below ("PATHS" section) !!



######## OUTPUTS :
### this script outputs many intermediary files (in the working repertory), you will mainly be interested in:
# 'orfs_contam': for each contaminant CDS : ID / corresponding scaffold / inferred taxonomic group (based on diamond blast results)
# 'tax_contam': for each contaminant scaffold : number of CDS / proportion of these contam CDS falling in each taxonomic group (-> detect chimeric scaffold, or problems in taxonomic assignment)
# 'HGT': ID of all HGT candidate CDS / inferred taxonomy of the donor species (based on diamond blast results) / species of the best diamond blast hit [! these candidates require a subsequent validation !]
# 'incertains': ID of all 'uncertain' CDS (i.e. suspicious, but alone on their scaffold -> synteny step not conclusive) / inferred taxonomy (based on diamond blast results) / species of the best diamond blast hit
# 'summ_after_gmap1': columns 1-14 correspond to intermediary results / col.15-34 correspond to final results; note that col.26-34 correspond to size distribution (min, max, average) of different types of scaffolds (contaminant / containing at least 1 potential HGT / uncertain)
# 'settings': file with the date of the analyses & the versions of the softwares.



######## REQUIREMENTS:
# requires a reference (custom) database used for the first diamond blast step. [you need to specify the absolute path; see below section "PATHS"]
# requires for each species / assembly: one file with the scaffolds and two files with the CDS (both RNA and protein sequences) -> 3 arguments required to run the script
# WARNING: headers should be consistent between the CDS_RNA and CDS_proteins fasta files !! (at least the first field of the header should be exactly similar in both files)
# required softwares: DIAMOND blast / Gmap / samtools



######## USAGE:
# Usage: foreign_cds_detection.sh [path/to/scaffold/file] [path/to/cds/file] [path/to/pep/file] [code_species] [number of CPU] [path/to/reference/database]
# this script should be launched from the folder in which all results will be written / output.
# ARG 1: fasta file with all scaffolds for this assembly
genome=$1
# ARG 2: fasta file with CDS sequences (RNA)
cds=$2
# ARG 3: fasta file with CDS sequences (proteins)
pep=$3
# ARG 4: a (short) code to be used to name all output files (species and/or assembly)
sp=$4
# ARG 5: number of CPU to be used for the DIAMOND blast and Gmap steps
cpu=$5
# ARG 6: absolute path to the reference database you want to use for the BLAST step (should be DIAMOND formatted, i.e. "database.dmnd"; omit the 'dmnd' extension):
customdb=$6



######## PATHS: !! TO BE ADAPTED TO YOUR ARCHITECTURE !!!!
# if DIAMOND blast and Gmap are not in your $PATH, you may need to specify their absolute path in the code below.
# to be noted:  this script copies the 3 fasta files (scaffold / cds / pep) in the working rep where the script was launched, and they are deleted at the end of the script (for cluster use).

######################################################################################################################################################################################################################



######## save the date of analyses and the software versions:
DATE=`date +%Y-%m-%d`
echo "Analyses were run the $DATE" > settings
echo `diamond --version` >> settings
echo `gmap --version | sed '3q;d'` >> settings


# function to retrieve all scaffolds linked to a CDS using the scaffold_bounded file (i.e. the scaffold on which this CDS mapped, as well as each scaffold linked to it)
get_all_scaffolds()
{

    orf=$1
    sc=`grep -w "$orf" "$sp"_gmap1/"$sp"_mapping | cut -f2 | head -n 1`

    if [[ -z $sc ]]; then
        echo ""
    else

        scaffs=`awk 'BEGIN{RS=">\n";FS="\n";} {for (i = 1; i<NF; i++){ if ($i == "'$sc'") {print $0 }; } }' scaffold_bounded_"$sp"`
        if [ "$scaffs" == "" ]; then
            echo "$sc"
        else
            echo "$scaffs"
        fi
    fi

}




############ START OF THE ACTUAL ANALYSES :


# first the files used by the pipeline are moved to the working directory (only useful for cluster use)
rep="$(pwd)"
cp "$genome" "$rep"/genome.fa
cp "$cds" "$rep"/cds.fa
cp "$pep" "$rep"/pep.fa	
genome="$(ls "$rep"/genome.fa)"
cds="$(ls "$rep"/cds.fa)"
pep="$(echo "$rep"/pep.fa)"

samtools faidx "$pep"
samtools faidx "$genome"

norf="$(grep -c '>' "$pep")"


######## 1. DIAMOND BLASTP step:

echo "Start the first DIAMOND BLAST step for $sp."

# first create the DIAMOND alignment archive (DAA) output (stores all information):
# option --masking : enable masking of low complexity regions (0/1=default)
diamond blastp -d $customdb -q "$pep" --threads $cpu --masking 1 --evalue 1E-5 --max-target-seqs 150 --more-sensitive -o "$sp"_diamond_blastp.daa --outfmt 100 1>"$sp"_diamond_blastp.log 2>"$sp"_diamond_blastp.err

# then create, from the DAA, the accurate tabular output: (here, we require 14 fields):
diamond view --daa "$sp"_diamond_blastp.daa --out "$sp"_diamond_blastp.tab --threads $cpu --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen sseq 1>"$sp"_diamond_view.log 2>"$sp"_diamond_view.err
# this outformat corresponds to the 12 default + 2 additional fields:
# qlen 		Query sequence length
# sseq 		Aligned part of subject sequence (here the gaps in the alignment are not indicated)

blast="$(echo "$sp"_diamond_blastp.tab)"

mkdir "$sp"_blast


# first filter: minimum of 40% identity, min 75AA alignment length and Evalue <= E-10 -> list of all CDS with at least 1 confident blast hit:
# sequential filter: first on Evalue (< 1e-10), then on %identity (>40), then on length(>75AA)
LC_ALL=C sort -k 11,11 -g $blast | LC_ALL=C awk '$11<=1e-10{print}' | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}'| cut -f1 | sort | uniq > "$sp"_blast/"$sp"_list_withhit
n_with_blast_hit="$(wc -l < "$sp"_blast/"$sp"_list_withhit)"

comm -3 <(grep -a '>' "$pep" | sed 's/>//g' | cut -d ' ' -f1 | sort) <(sort "$sp"_blast/"$sp"_list_withhit) > "$sp"_blast/"$sp"_nohit
n_no_blast_hit="$(wc -l < "$sp"_blast/"$sp"_nohit)"


# prepare the files for each taxonomic group (see 'categ_all') and remove them if they already exist:
#  2>/dev/null redirects the error message (if the file does not exist) to the 'dump' folder = avoid polluting the screen/log.
categ_all="archaea eubacteria fungi plant protist uncertainonlynm arthropod othermetazoa"
for categ in $categ_all
do
	rm "$sp"_blast/"$sp"_"$categ" 2>/dev/null; touch "$sp"_blast/"$sp"_"$categ"
done


####### parse the BLAST output, based on defined criteria for taxonomic assignment:
# on the 10 best unique hits (unique -> keep only 1 hit per reference sequence, as to not give too much weight to a single seq which might be erroneous).
# hits from at least 2 species are required (otherwise: discard this CDS)
# we choose to be stringent for the 'arthropod' genes, because later used to confirm that 'dubious' genes are actually integrated IN the host genome]
# arthropod = only metazoa hits and at least 70% arthropod (including at least 2 arthropod species)
# [DEPRECATED] othermetazoa (actually: this set is not used afterwards) = at least 70% of 'othermetazoa' hits and NO arthropod hit ! [DEPRECATED] 
# for other groups: criteria are relaxed to account for potential contamination (from symbionts, diet, or other) in genome assemblies.
# plant / eubacteria / archaea / fungi / protist = at least 70% of the corresponding taxonomic group
# uncertain = all CDS meeting none of the above requirements for a clear taxonomic assignment
# (including) uncertainonlynm = same as above BUT only non-metazoa hits.


# loop on all CDS with at least 1 confident blast hit ( with 1e-10 evalue, 40% identity and alignment length of 75 AA at least):

for contig in `cat "$sp"_blast/"$sp"_list_withhit`
do
	
	# 1. create 'contigtmp' with the 10 best unique hits that we will be using:
	# awk '!x[$0]++' enables to remove the redundant ref sequences while keeping the order (by evalue) = keeping the best hit of this ref sequence
	grep $contig $blast | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}'| LC_ALL=C sort -k 11,11 -g | LC_ALL=C awk '$11<=1e-10{print}' | cut -f2 | awk '!x[$0]++' | head -n 10 > contigtmp
	# list of the 'detailed group' of the 10 hits (ie arthropod / eubacteria etc)
	cut -d '|' -f2 contigtmp > temp #->categ
	# uniq ( metazoa / non-metazoa ) of the 10 hits:
	cut -d '|' -f1 contigtmp | sort | uniq > temp2
	# cut the 'detailed | species' field  (ex: arthropod|drosophila_melanogaster)
	cut -d '|' -f2,3 contigtmp | sort | uniq > temp3
	
	occ=`cat temp | wc -l`
	fifty=`awk "BEGIN {print "$occ"/2}"`
	seventy=`awk "BEGIN {print "$occ"*7/10}"`
	
	# set all 'counters' at 0 [ in the 10 best hits: arthropod = nber of arthropod HITS / narthropod= number of arthropod SPECIES ]
		arthropod="$(echo "0")"; narthropod="$(echo "0")"
	
	# number of different species in the 10 best hits:
	nspec="$(cut -d '|' -f2 temp3 | sort | uniq | wc -l)"
	
	# CASE: arthropod + othermetazoa
	if  grep -Fxq "arthropod" temp; then
	    arthropod=`sort temp | uniq -c | grep 'arthropod' | awk '{print $1}'`
	fi
	
	narthropod="$(grep -w 'arthropod' temp3 | cut -d '|' -f2 | sort | uniq | wc -l)"
	
	if  [[ `uniq temp` = "arthropod" ]] && (( $(echo "$narthropod >= 2" | bc -l) )); then
	    echo "$contig" >>"$sp"_blast/"$sp"_arthropod
	fi
	if [[ `cat temp2` = "metazoa" ]] && (( $(echo "$arthropod >= $seventy" | bc -l) )) && grep -Fxq "othermetazoa" temp; then
	    echo "$contig" >> "$sp"_blast/"$sp"_arthropod
	fi
	
	
	if  grep -Fxq "othermetazoa" temp; then
	    doubt=`sort temp | uniq -c | grep -w "othermetazoa" | awk '{print $1}'`
	    ndoubt="$(grep -w "othermetazoa" temp3 | cut -d '|' -f2 | sort | uniq | wc -l)"
	    if [[ `uniq temp` = "othermetazoa" ]] && (( $(echo "$ndoubt >= 2" | bc -l) )) && [[ "$narthropod" -eq 0 ]]; then
			echo "$contig" >> "$sp"_blast/"$sp"_othermetazoa
	    fi
	    if (( $(echo "$doubt >= $seventy" | bc -l) )) && (( $(echo "$ndoubt >= 2" | bc -l) )) && [[ `uniq temp` != "othermetazoa" ]] && [[ "$narthropod" -eq 0 ]]; then
			echo "$contig" >> "$sp"_blast/"$sp"_othermetazoa
	    fi
	fi
	
	
	## CASE = suspected FOREIGN : loop on the different 'categ_pure' of FOREIGN ( = archaea / eubacteria / fungi / plant / protist ) :
	categ_pure="archaea eubacteria fungi plant protist"
	for categ in $categ_pure
	do
	    # set the counters at 0 for each category:
	    doubt="$(echo "0")"; ndoubt="$(echo "0")"
	
	    if  grep -Fxq "$categ" temp; then
		
			doubt=`sort temp | uniq -c | grep -w "$categ" | awk '{print $1}'`
			ndoubt="$(grep -w "$categ" temp3 | cut -d '|' -f2 | sort | uniq | wc -l)"
			if  [[ `uniq temp` = "$categ" ]] && (( $(echo "$ndoubt >= 2" | bc -l) )); then
			    echo "$contig" >>"$sp"_blast/"$sp"_"$categ"
			fi
			if (( $(echo "$doubt >= $seventy" | bc -l) )) && (( $(echo "$ndoubt >= 2" | bc -l) )) && [[ `uniq temp` != "$categ" ]]; then
			    echo "$contig" >> "$sp"_blast/"$sp"_"$categ"
			fi
	
	    fi
	
	done

done


# create the list of CDS with no 'reliable' taxonomic assignment = "uncertain"
cat "$sp"_blast/"$sp"_arthropod "$sp"_blast/"$sp"_othermetazoa > "$sp"_blast/"$sp"_all_tax_assigned
for categ in $categ_pure
do
	cat "$sp"_blast/"$sp"_"$categ" >> "$sp"_blast/"$sp"_all_tax_assigned
done

comm -3 <(sort "$sp"_blast/"$sp"_list_withhit) <(sort "$sp"_blast/"$sp"_all_tax_assigned) > "$sp"_blast/"$sp"_uncertain


# then, create the sub-list of these 'unreliable-taxonomy' CDS which have ONLY non-metazoa hits = "uncertainonlynm" => will be kept
for contig in `cat "$sp"_blast/"$sp"_uncertain`
do
	nspec="$(echo "0")"
	grep $contig $blast | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}'| LC_ALL=C sort -k 11,11 -g | LC_ALL=C awk '$11<=1e-10{print}'| cut -f2 | awk '!x[$0]++' | head -n 10 > contigtmp
	cut -d '|' -f1 contigtmp | sort | uniq > temp2
	cut -d '|' -f2,3 contigtmp | sort | uniq > temp3
	nspec="$(cut -d '|' -f2 temp3 | sort | uniq | wc -l)"
	if [[ `cat temp2` = "nonmetazoa" ]] && (( $(echo "$nspec >= 2" | bc -l) )); then
	    echo "$contig" >> "$sp"_blast/"$sp"_uncertainonlynm
	fi
done

# cleaning the temp files:
rm contigtmp 2>/dev/null; rm temp 2>/dev/null; rm temp2 2>/dev/null; rm temp3 2>/dev/null


# also generate a summary file "$sp"_all_tax with the taxonomic assignment for all CDS (when it was achievable):
rm "$sp"_all_tax 2>/dev/null; touch "$sp"_all_tax
# for the taxonomic groups: archaea / eubacteria / fungi / plant / protist / uncertainonlynm /arthropod/ othermetazoa
for categ in $categ_all
do
	for orf in `cat "$sp"_blast/"$sp"_"$categ"`
	do
	    printf "$orf\t$categ\n" >> "$sp"_all_tax
	done
done
# then, remove all 'uncertain', because tricky to manipulate afterwards ! (doublons with 'uncertainonlynm')
sed -i '/uncertain$/d' "$sp"_all_tax


# summarize these intermediary results (these will be printed afterwards in the joint summary file 'summ_after_gmap1' with blast & Gmap intermediary results.)
# WARNING: these figures should not be used as final results, but could be used for debugging)
norfwithtax="$(wc -l < "$sp"_blast/"$sp"_all_tax_assigned)"
norfarthropod="$(wc -l < "$sp"_blast/"$sp"_arthropod)"
o="$(wc -l < "$sp"_blast/"$sp"_othermetazoa)"
norfmetazoa=`awk "BEGIN {print "$norfarthropod"+"$o"}"`
nhgteub="$(wc -l < "$sp"_blast/"$sp"_eubacteria)"
nhgtarch="$(wc -l < "$sp"_blast/"$sp"_archaea)"
nhgtprot="$(wc -l < "$sp"_blast/"$sp"_protist)"
nhgtfung="$(wc -l < "$sp"_blast/"$sp"_fungi)"
nhgtplant="$(wc -l < "$sp"_blast/"$sp"_plant)"
nhgt=`awk "BEGIN {print "$nhgteub"+"$nhgtarch"+"$nhgtprot"+"$nhgtfung"+"$nhgtplant"+"$o"}"`
nhgtunconlynm="$(wc -l < "$sp"_blast/"$sp"_uncertainonlynm)"


echo "End of the BLAST step for $sp."


################ 2. Gmap step (synteny):

echo "Start the Gmap step for $sp."

mkdir "$sp"_gmap1


#### 2a. for Gmap, first need to 'index' the genome assembly (here, one genomic scaffold per FASTA entry):
# redirect the standard output / error to "$sp"_refgenome/gmap_build.log and gmap_build.err
echo "Start indexing the genome assembly $sp with/for Gmap"
cd "$rep"; mkdir "$sp"_refgenome; cd "$sp"_refgenome
gmap_build -D "$rep"/"$sp"_refgenome -d "$sp" "$genome" 1>gmap_build.log 2>gmap_build.err

cd "$rep"



#### 2b. map all CDS (for which taxonomic assignment was achievable / reliable in the previous blast step) to the ref genomic scaffolds:
# we have defined a minimum length of 100 bp and min identity of 95% for this alignment step

# create the list of to-be-Gmapped CDS (for which taxonomic assignment was achievable in the previous blast step), if it doesn't exist yet:
cat "$sp"_blast/"$sp"_all_tax_assigned "$sp"_blast/"$sp"_uncertainonlynm > "$sp"_gmap1/"$sp"_orf_with_tax


# retrieve all the corresponding cDNA sequences (=> 'all_tax_assigned' + 'uncertainonlynm') in a fasta file:
cat "$sp"_gmap1/"$sp"_orf_with_tax | xargs -n 1 samtools faidx "$cds" > "$sp"_cds


# Gmap: map all cDNA on the 'indexed' genomic scaffolds (min 95% identity).
# default: maximum 5 output paths for each mapped cDNA (multiple paths + chimeric alignments = ESTs whose 5 and 3 ends map to different genomic regions.)
# so here, option --npaths=0 to get ONLY a single alignment plus chimeric alignments, for each cDNA (exclude multiple alignements here)
# we keep chimeric alignments as could be due to the fragmentation of genome assembly = convey some information.
# We redirect the SAM std output to "$sp"_gmap.sam and std error to "$sp"_gmap.err
gmap -D "$rep"/"$sp"_refgenome -d "$sp" -B 5 -t $cpu -f samse --no-sam-headers --npaths=0 --min-identity=0.95 "$sp"_cds 1>"$sp"_gmap1/"$sp"_gmap.sam 2>"$sp"_gmap1/"$sp"_gmap.err

# cleaning:
rm "$sp"_cds

# Reminder: to check / visualize the alignment (in SAM format):use sam2pairwise to see in a human-readable view the alignment (open the output with geany)
#grep 'Tsi_TRINITY_DN29652_c0_g3_i2|m.6970' ex_sam | ~/software/sam2pairwise-master/src/sam2pairwise > ~/tmp; geany ~/tmp



#### 2c. create a table summarizing the info for each mapped cDNA:  CDS id / scaffold on which it mapped / taxonomy of the CDS
# when no mapping was achievable for a given CDS: indicated in the SAM file by $2==4 AND $3 = '*'
# require a min. length of 100 bp (of the aligned segments (~ exons)) (~ exon length in metazoa & empirical on some visual tests)

rm "$sp"_gmap1/"$sp"_nomapping 2>/dev/null; touch "$sp"_gmap1/"$sp"_nomapping
rm "$sp"_gmap1/"$sp"_mapping 2>/dev/null; touch "$sp"_gmap1/"$sp"_mapping


for orf in `cat "$sp"_gmap1/"$sp"_orf_with_tax`
do	
    tax="$(echo NA)"
    grep -w "$orf" "$sp"_gmap1/"$sp"_gmap.sam > samtemp
    # if the CDS was not reliably mapped:
    if  [[ `cut -f2 samtemp` == 4 ]]; then
		echo "$orf" >> "$sp"_gmap1/"$sp"_nomapping
    else
		# assign the taxonomic group to each CDS based on the blast step ("$sp"_all_tax)
		# here the taxonomic groups are: archaea / eubacteria / fungi / plant / protist / arthropod / othermetazoa / uncertainonlynm
		tax="$(grep -w "$orf" "$sp"_all_tax | cut -f2)"
	
	
		## now we summarize this mapping/tax info (for each CDS), in the file "$sp"_mapping:
		# here, applies to the cases where the CDS mapped on 1 or more scaffolds:
		nmap="$(cat samtemp | wc -l)"
		for i in $(seq 1 1 $nmap)
		do
		    # retrieve the focal line:
		    sed "${i}q;d" samtemp > subtemp
		    # then, retrieve the MD string to calculate the number of aligned bases (matches=md1 + mismatches=md2):
		    # attention! the SAM output of Gmap is weird: when a given cDNA maps to > 1 scaffold, the MD tag is in the 13th field instead of the 12th...
		    md="$(tr '[ \t]' '[\n\n]' < subtemp | grep '^MD' | cut -d':' -f3)"
		    # md1 = matches = sum of the numbers indicative of matches in the MD string
		    md1="$(echo $md | grep -o -E '[0-9]+'| awk '{ SUM += $1} END { print SUM }')"
		    # md2 = mismatches = sum of the letters indicative of mismatches in the MD string
		    md2="$(echo $md | sed 's/[0-9]//g' | wc -c)"
		    # require a min length of alignment (for the mapping parts ~ exons) of 100 bp:
		    if (( $(echo "$md1 + $md2 - 1 >= 100" | bc -l) )); then
				scaff="$(cut -f3 subtemp)"
				# print one line per scaffold on which this CDS mapped: CDS id / scaffold ID / taxonomy assigned to this CDS
				printf "$orf\t$scaff\t$tax\n" >> "$sp"_gmap1/"$sp"_mapping
		    fi
		    rm subtemp
		done
    fi
    rm samtemp
done



# loop to check if there is any 'NA' taxonomy remaining => problem
if grep -Fwxq "NA" <(cut -f3 "$sp"_gmap1/"$sp"_mapping); then
    echo "!! There is a problem of taxonomic assignment for some CDS in the species $sp [Gmap step] !!"
fi



###### scaffolding using the CDS mapping/linking several scaffolds (=  called 'chimeric' in Gmap) to retrieve more HGT candidates (= convey the 'arthropod' tag to linked scaffolds):
# if a cDNA maps on several scaffolds, then they could be reunited in the assembly. a file is created, scaffold_bounded_specie with every group of scaffolds
# that could be reunited
# (this is an iterative process (scaff can be linked by several cDNA), that's why the loop is so complicated)

rm scaffold_bounded_"$sp" 2>/dev/null; touch scaffold_bounded_"$sp"

liste_orfs_liant=`cat "$sp"_gmap1/"$sp"_mapping | cut -f1 | sort | uniq -d `

until [[ -z $liste_orfs_liant ]] 
do
    tab_orfs=($liste_orfs_liant)
    orf=${tab_orfs[0]}

    scaff_orf=`grep -w "$orf" "$sp"_gmap1/"$sp"_mapping | cut -f2 | sort | uniq `

    tout_orfs_lies=`grep -w "$scaff_orf" "$sp"_gmap1/"$sp"_mapping | cut -f1  | sort | uniq `
    orfs_lies=`comm -1 -2 <( echo "$tout_orfs_lies") <(echo "$liste_orfs_liant") | sort `	    
    group_scaff=""
    
    until [[ -z $orfs_lies ||  -z $liste_orfs_liant ]]
    do		
		scaff_lies=`grep "$orfs_lies" "$sp"_gmap1/"$sp"_mapping | cut -f2 | sort | uniq `
		group_scaff=`printf "$group_scaff\n$scaff_lies"`
		liste_orfs_liant=`comm -2 -3 <(echo "$liste_orfs_liant") <(echo "$orfs_lies") | sort | uniq `
		tout_orfs_lies=`grep -w "$scaff_lies" "$sp"_gmap1/"$sp"_mapping | cut -f1 | sort | uniq `
		orfs_lies=`comm -1 -2 <(echo "$tout_orfs_lies") <(echo "$liste_orfs_liant") | sort `
    done
    
    group_scaff=`echo "$group_scaff" | sort | uniq`	    
    printf ">$group_scaff\n" >> scaffold_bounded_"$sp"
done




############ now, we can assign each foreign candidate (detected in the 1st blast step) to one of these groups: contam / potential HGT / uncertain :

#### we start by extracting HGT candidates from all foreign candidates (based on synteny):

# keep only HGT candidates = those which are located on (or linked to) 'arthropod-like' scaffolds (= on which at least 1 'arthropod-like' CDS mapped)
categ_HGT2="archaea eubacteria fungi plant protist uncertainonlynm"
# create the list of foreign candidates (after the blast step, including uncertainonlynm), if it doesn't exist yet:
rm "$sp"_gmap1/"$sp"_candid_after_blast2 2>/dev/null; touch "$sp"_gmap1/"$sp"_candid_after_blast2
for cat in $categ_HGT2
do
    cat "$sp"_blast/"$sp"_"$cat" >> "$sp"_gmap1/"$sp"_candid_after_blast2
done

# prepare the files for each category (archaea / eubacteria / fungi / plant / protist / uncertainonlynm) :
for categ in $categ_HGT2
do
    rm "$sp"_gmap1/"$sp"_"$categ" 2>/dev/null; touch "$sp"_gmap1/"$sp"_"$categ"
done
		
# create a new version of sp_mapping without any 'othermetazoa' [reminder: 'othermetazoa' is a deprecated category]:
grep -vw 'othermetazoa' "$sp"_gmap1/"$sp"_mapping > "$sp"_gmap1/"$sp"_mapping2	

### scann all candidates, decide if they are to be kept as HGT or not, and assign them to their taxonomic group (eubacteria / plant / ...):
for cand in `cat "$sp"_gmap1/"$sp"_candid_after_blast2`
do
    # retrieve the taxonomic group of this HGT-candidate:
    taxcand="$(grep -w "$cand" "$sp"_all_tax | cut -f2 )"
    # on which scaffold(s) did it map or is it linked (can be more than 1): (function get_all_scaffolds() defined on top of this script)
    scaff2="$(get_all_scaffolds $cand)"
    # then retrieve the taxonomic assignment of ALL ORFs which have been mapped to this (or these) given scaffold(s):
    cat "$sp"_gmap1/"$sp"_mapping2 | grep -w "$scaff2" | cut -f3 | sort | uniq > temp
    ## now, we extract the candidates which shared at least 1 scaffold with another 'arthropod-like' gene => HGT candidates [ !! contamination filter, sensitive to assembly N50 !! ]
    if (grep -Fwxq "arthropod" temp); then
		echo "$cand" >> "$sp"_gmap1/"$sp"_"$taxcand"
    fi
    rm temp
done


# create the list of all HGT-candidates (including 'uncertainonlynm') surviving this Gmap filter: "$sp"_candid_after_gmap2:
categ_HGT2="archaea eubacteria fungi plant protist uncertainonlynm"
rm "$sp"_gmap1/"$sp"_candid_after_gmap2 2>/dev/null; touch "$sp"_gmap1/"$sp"_candid_after_gmap2
for ct in $categ_HGT2
do
    cat "$sp"_gmap1/"$sp"_"$ct" >> "$sp"_gmap1/"$sp"_candid_after_gmap2
done

# summarize these results to print them in the summary file (later):
n_orf_mapped="$(cut -f1 "$sp"_gmap1/"$sp"_mapping | sort | uniq | wc -l)"
n_orf_to_be_mapped="$(wc -l < "$sp"_gmap1/"$sp"_orf_with_tax)"
p_orf_mapped="$(echo $(( n_orf_mapped *100 / n_orf_to_be_mapped )))"
nhgteub2="$(wc -l < "$sp"_gmap1/"$sp"_eubacteria)"
nhgtarch2="$(wc -l < "$sp"_gmap1/"$sp"_archaea)"
nhgtprot2="$(wc -l < "$sp"_gmap1/"$sp"_protist)"
nhgtfung2="$(wc -l < "$sp"_gmap1/"$sp"_fungi)"
nhgtplant2="$(wc -l < "$sp"_gmap1/"$sp"_plant)"
nhgt2=`awk "BEGIN {print "$nhgteub2"+"$nhgtarch2"+"$nhgtprot2"+"$nhgtfung2"+"$nhgtplant2"}"`
nhgtunconlynm2="$(wc -l < "$sp"_gmap1/"$sp"_uncertainonlynm)"




####### then, extract all contam CDS, and generate some summary / stats files :


# putative taxonomy for HGT candidates (general taxonomy as estimated by the diamond blastp step):
printf "ID_orf\ttax_orf\tbest_hit\n" > HGT_"$sp"

# putative taxonomy for uncertain genes = CDS alone on their scaffold (general taxonomy as estimated by the diamond blastp step):
printf "ID_orf\ttax_orf\tbest_hit\n" > incertains_"$sp"
	
# summary file about the taxonomy estimation (by blast) of all CDS in each contaminant scaffold :
#-> !! If there are scaffold mixing different taxonomy (eg 50% CDS that are attributed to plant and 50% that are eubacteria) the genome assembly might not be reliable (eg chimeras)
printf "scaffold\tnb_orfs\teubacteria\tprotist\tarchaea\tfungi\tplant\tonlynm\n" > tax_contam_"$sp"

# (intermediary / working) file listing all the contaminant CDS, their taxonomy (blast) and the scaffold to which they map
rm orfs_contam_"$sp" 2>/dev/null; touch orfs_contam_"$sp"	

# variables to run stats on the 'foreign' scaffolds (will be written in the summ_after_gmap1 file)
dna_fai="$(ls "$genome".fai)"

blast="$(echo "$sp"_diamond_blastp.tab)"
info="$(echo "$sp"_gmap1)"
min_tscaff_contam=30000000000000
med_tscaff_contam=0
max_tscaff_contam=0
min_tscaff_uncertain=30000000000000
med_tscaff_uncertain=0
max_tscaff_uncertain=0	
min_tscaff_hgt=30000000000000
med_tscaff_hgt=0
max_tscaff_hgt=0
nb_scaff_hgt=0

candid_blast=`sort "$info"/"$sp"_candid_after_blast2`
HGT=`sort "$info"/"$sp"_candid_after_gmap2`
candid_contamination=`comm -2 -3 <( echo "$candid_blast") <(echo "$HGT") | sort`

if [ -z "$candid_contamination" ]; then
    nb_tot_nhgt=0
else
    nb_tot_nhgt=`cat <(echo "$candid_contamination") | wc -l`	    
fi

conta_mapping=`grep -w "$candid_contamination" "$info"/"$sp"_mapping2 | cut -f1 | sort | uniq`

if [ -z "$conta_mapping" ]; then
    nb_mapping=0
else	    
    nb_mapping=`cat <(echo "$conta_mapping") | wc -l`
fi

nb_nomapping=$((nb_tot_nhgt - nb_mapping))	# correspond to the 'suspicious' candidates which did not map (reliably) on the genome assembly.


# analyses on the scaffolds	carrying the HGT candidates (size distribution):
scafs="$(grep -w "$HGT" "$info"/"$sp"_mapping2 |cut -f2|sort|uniq)"

for scaffold in $scafs 
do		    
    taille_scaffold="$(echo `grep -w "$scaffold" "$dna_fai" | cut -f2`)"
   	    
    if [[ "$taille_scaffold" -lt "$min_tscaff_hgt" ]] ; then		
		min_tscaff_hgt="$taille_scaffold"		
    fi
    
    if [[ "$taille_scaffold" -gt "$max_tscaff_hgt" ]] ; then		
		max_tscaff_hgt="$taille_scaffold"
    fi
    
    med_tscaff_hgt=$((med_tscaff_hgt + taille_scaffold))
    nb_scaff_hgt=$(($nb_scaff_hgt + 1))  
done

if [ $nb_scaff_hgt -gt 1 ]; then
    med_tscaff_hgt=$((med_tscaff_hgt / nb_scaff_hgt))
fi



# analyses on the HGT candidates at the CDS level:	
for orf in $HGT
do	    
    tax=`grep -w "$orf" "$info"/"$sp"_mapping | cut -f3 | sort -u`
    best_hit=`grep -w $orf $blast | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}'| LC_ALL=C sort -k 11,11 -g | LC_ALL=C awk '$11<=1e-10{print}' | cut -f2 | awk '!x[$0]++' | head -n 1 `
    specie_field=`echo "$best_hit" | cut -d "|" -f3`    
    printf "$orf\t$tax\t$specie_field\n" >> HGT_"$sp"  	   	         
done




# analyses on the scaffolds with only 1 suspicious CDS (=uncertain) and those with >1 (and only) non-arthropod CDS (=contaminations)
scaffolds_candid_conta=`grep -w "$conta_mapping" "$info"/"$sp"_mapping2 | cut -f2 | sort | uniq`

if [ -z "$scaffolds_candid_conta" ]; then
    nb_scaff_nhgt=0
else
    nb_scaff_nhgt=`wc -l <( echo "$scaffolds_candid_conta") | cut -d ' ' -f1`
fi

for sc in $scaffolds_candid_conta
do
    taille_scaffold="$(echo `grep -w "$sc" "$dna_fai" | cut -f2`)"   
    nb_orfs=`grep -c -w "$sc" "$info"/"$sp"_mapping2 | sort -u`
    
    # analyses on the suspicious CDS alone on the scaffold = "uncertain"
    if [ $nb_orfs -eq 1 ]; then	
	    if [[ "$taille_scaffold" -lt "$min_tscaff_uncertain" ]] ; then 
			min_tscaff_uncertain="$taille_scaffold"		
	    fi	
	    if [[ "$taille_scaffold" -gt "$max_tscaff_uncertain" ]] ; then		
			max_tscaff_uncertain="$taille_scaffold"
	    fi  
	    med_tscaff_uncertain=$((med_tscaff_uncertain + taille_scaffold))   			
		orf_line=`grep -w "$sc" "$info"/"$sp"_mapping2 | sort -u`
		orf=` echo "$orf_line" | cut -f1`
		tax_orf=` echo "$orf_line" | cut -f3`
		best_hit=`grep -w $orf $blast | LC_ALL=C sort -k 3,3 -g | LC_ALL=C awk '$3>=40{print}' | LC_ALL=C sort -k 4,4 -g | LC_ALL=C awk '$4>=75{print}'| LC_ALL=C sort -k 11,11 -g | LC_ALL=C awk '$11<=1e-10{print}' | cut -f2 | awk '!x[$0]++' | head -n 1 `
		specie_field=`echo "$best_hit" | cut -d "|" -f3 `  
	    printf "$orf\t$tax_orf\t$specie_field\n" >> incertains_"$sp"	    	    			       		
	
    else
		# analyses on the scaffolds with >1 CDS non-arthropod (=contaminations)
		if [[ "$taille_scaffold" -lt "$min_tscaff_contam" ]] ; then 
			min_tscaff_contam="$taille_scaffold"		
	    fi	
	    if [[ "$taille_scaffold" -gt "$max_tscaff_contam" ]] ; then		
			max_tscaff_contam="$taille_scaffold"
	    fi  
	    med_tscaff_contam=$((med_tscaff_contam + taille_scaffold)) 
			
		# variables used for the taxonomy summary file that sums up the taxonomy of the contaminant scaffolds.
		tax_eubact=0
		tax_prot=0
		tax_arch=0
		tax_fun=0
		tax_pl=0
		tax_onlynm=0	
		orfs=`grep -w "$sc" "$info"/"$sp"_mapping2 | cut -f1 | sort -u`
		nb=0
			
		for orf in $orfs
		do	
			tax=`grep -w "$orf" "$info"/"$sp"_mapping2 | cut -f3 | sort -u `
			printf "$orf\t$sc\t$tax\n" >> orfs_contam_"$sp"
			case "$tax" in
				"eubacteria") tax_eubact=$(($tax_eubact + 1));;
				"protist") tax_prot=$(($tax_prot + 1));;
				"archaea") tax_arch=$(($tax_arch + 1));;
				"fungi") tax_fun=$(($tax_fun + 1));;
				"plant") tax_pl=$(($tax_pl +1));;
				*) tax_onlynm=$(($tax_onlynm + 1));;			 
		    esac
			    nb=$(($nb + 1))
					    	   
		done
		# then, fill the (more general) table summarizing the taxonomy of the CDS on all contaminant scaffolds:
		tax_eubact=$(($tax_eubact * 100))
		tax_eubact=$((tax_eubact / nb))
		tax_prot=$(($tax_prot * 100))
		tax_prot=$((tax_prot / nb))
		tax_arch=$(($tax_arch * 100))
		tax_arch=$((tax_arch / nb))		
		tax_fun=$(($tax_fun * 100))
		tax_fun=$((tax_fun / nb))
		tax_pl=$(($tax_pl * 100))
		tax_pl=$((tax_pl / nb))
		tax_onlynm=$(($tax_onlynm * 100))
		tax_onlynm=$((tax_onlynm / nb))	
		printf "$sc\t$nb\t$tax_eubact\t$tax_prot\t$tax_arch\t$tax_fun\t$tax_pl\t$tax_onlynm\n" >> tax_contam_"$sp"
		#REMINDER: printf "scaffold\tnb_orfs\teubacteria\tprotist\tarchaea\tfungi\tplant\tonlynm\n" > tax_contam_"$sp"
    fi
        
done


nb_incertains=`cat incertains_"$sp" | wc -l`
nb_incertains=$(($nb_incertains - 1))
nb_contam=$((nb_mapping - nb_incertains))
nb_scaff_contam=$((nb_scaff_nhgt - nb_incertains))
	
if [ $nb_scaff_contam -gt 1 ]; then
    med_tscaff_contam=$((med_tscaff_contam / nb_scaff_contam))
fi
if [ $nb_incertains -gt 1 ]; then
    med_tscaff_uncertain=$((med_tscaff_uncertain / nb_incertains))
fi
if [[ $min_tscaff_hgt -eq 30000000000000 ]]; then
    min_tscaff_hgt=0
fi
if [[ $min_tscaff_contam -eq 30000000000000 ]]; then
    min_tscaff_contam=0
fi
if [[ $min_tscaff_uncertain -eq 30000000000000 ]]; then
    min_tscaff_uncertain=0
fi	


# print all results (intermediary col1-14 / final col15-34)
printf ""species"\t"N_orf"\t"N_orf_with_tax"\t"N_hgt1"\t"N_hgt1_eubact"\t"N_hgt1_arch"\t"N_hgt1_protist"\t"N_hgt1_fungi"\t"N_hgt1_plant"\t"N_hgt1_uncertainonlynm"\t"N_orf_arthropod"\t"N_orf_metazoa"\t"N_orf_mapped"\t"percent_orf_mapped"\t"N_hgt2"\t"N_hgt2_eubact"\t"N_hgt2_arch"\t"N_hgt2_protist"\t"N_hgt2_fungi"\t"N_hgt2_plant"\t"N_hgt2_uncertainonlynm"\t"N_orf_nomapping"\t"N_orf_contamination"\t"N_orf_incertains"\t"N_scaff_contam"\t"min_tscaff_contam"\t"avg_tscaff_contam"\t"max_tscaff_contam"\t"min_tscaff_uncertain"\t"avg_tscaff_uncertain"\t"max_tscaff_uncertain"\t"min_tscaff_hgt"\t"avg_tscaff_hgt"\t"max_tscaff_hgt"\n" > "$rep"/summ_after_gmap1
printf "$sp\t$norf\t$norfwithtax\t$nhgt\t$nhgteub\t$nhgtarch\t$nhgtprot\t$nhgtfung\t$nhgtplant\t$nhgtunconlynm\t$norfarthropod\t$norfmetazoa\t$n_orf_mapped\t$p_orf_mapped\t$nhgt2\t$nhgteub2\t$nhgtarch2\t$nhgtprot2\t$nhgtfung2\t$nhgtplant2\t$nhgtunconlynm2\t$nb_nomapping\t$nb_contam\t$nb_incertains\t$nb_scaff_contam\t$min_tscaff_contam\t$med_tscaff_contam\t$max_tscaff_contam\t$min_tscaff_uncertain\t$med_tscaff_uncertain\t$max_tscaff_uncertain\t$min_tscaff_hgt\t$med_tscaff_hgt\t$max_tscaff_hgt\n" >> "$rep"/summ_after_gmap1


echo "End of the Gmap step for $sp"
cd "$rep"
rm pep.fa; rm pep.fa.fai  #=> for cluster use (to avoid big amount of data in current work directory)
rm genome.fa; rm genome.fa.fai
rm cds.fa; rm cds.fa.fai
rm "$sp"_gmap1/"$sp"_mapping2

## End of the analyses

