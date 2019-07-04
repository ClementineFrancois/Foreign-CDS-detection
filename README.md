# Foreign-CDS-detection

# About

This script identifies contaminant CDS / scaffolds in a genome assembly, without any a priori on the source(s) of contamination. It also outputs a list of potential HGT candidates, which require a subsequent (phylogenetic) validation. This script focuses on non-metazoan sources of contaminants and HGT, classified in five taxonomic groups: eubacteria, archaea, fungi, viridiplantae and ‘protists’. It was optimized and benchmarked in arthropods, but can be applied to any other taxa, provided that an adequate reference database is available.  

The first step of the script consists in a preliminary taxonomic assignment of each CDS based on sequence similarity (DIAMOND BLASTP against a protein reference database). CDS assigned to a non-metazoan group (eubacteria; archaea; viridiplantae; fungi; protists) are considered as “foreign” candidates. In addition, CDS assigned to Arthropoda are labelled as "confident-arthropod" to be later used in the script. The 10 best blast hits are considered instead of just the best one, as a way to account for potential contaminations and other sources of taxonomic mis-assignment in the reference database.

The second step of the script is a test of synteny. All foreign CDS candidates as well as the “confident-arthropod” CDS are mapped onto the genomic scaffolds using GMAP. A candidate is considered as a potential HGT if it was physically linked to (i.e., mapped to the same scaffold as) at least one “confident-arthropod” CDS. A candidate is considered as a contaminant if it mapped to a scaffold to which no “arthropod-confident” CDS mapped, and at least another non-metazoa CDS mapped. A candidate is considered as “uncertain” if it did not reliably map to any scaffold or if it was the only CDS to map to a given scaffold.

As an example, the genome assembly of Aedes aegypti (EnsemblMetazoa) was processed in 3.5 hours using 70 CPU.



# Data

- a genome assembly (scaffolds)

- the corresponding set of predicted coding sequences (CDS) (both as cDNA and protein sequences)
WARNING: headers should be consistent between the CDS_cDNA and CDS_protein fasta files (at least the first field of the header should be exactly similar in both files).

- a protein reference database (for the DIAMOND blast step). This database needs to be balanced and to represent the actual diversity of living organisms.

- required softwares: DIAMOND blast / Gmap / samtools


# Usage

foreign_cds_detection.sh [path/to/scaffold/file] [path/to/cds/file] [path/to/pep/file] [code_species] [number of CPU] [path/to/reference/database]

Example: ```foreign_cds_detection.sh assembly.fa cds.fa pep.fa aedes 50 /home/user/refdatabase```

This script should be launched from the folder in which all results will be written. 

Requires 6 arguments:
- ARG 1: path to the scaffolds file
- ARG 2: path to the fasta file with CDS sequences (cDNA)
- ARG 3: path to the fasta file with CDS sequences (protein)
- ARG 4: a (short) code to be used to name all output files
- ARG 5: number of CPU to be used for the DIAMOND blast and Gmap steps
- ARG 6: path to the reference database for the BLAST step (should be DIAMOND-formatted, i.e. "database.dmnd"; omit the 'dmnd' extension)



# To be noted

- This script copies the 3 fasta files (scaffold / cds / pep) in the working rep where the script was launched, and they are deleted at the end of the script (for cluster use).

- If DIAMOND blast and Gmap are not in your $PATH, you need to specify their absolute path in the code.



# Outputs

This script outputs many intermediary files (in the working repertory), you will mainly be interested in:

- 'orfs_contam': for each contaminant CDS : ID / corresponding scaffold / inferred taxonomic group (based on diamond blast results)

- 'tax_contam': for each contaminant scaffold : number of CDS / proportion of these contaminant CDS falling in each taxonomic group (-> to detect chimeric scaffolds, or problems in taxonomic assignment)

- 'HGT': ID of all HGT candidate CDS / inferred taxonomy of the donor species (based on diamond blast results) / species of the best diamond blast hit [! these candidates require a subsequent validation !]

- 'incertains': ID of all 'uncertain' CDS (i.e. suspicious, but alone on their scaffold -> synteny step was not conclusive) / inferred taxonomy (based on diamond blast results) / species of the best diamond blast hit

- 'summ_after_gmap1': columns 1-14 correspond to intermediary results / col.15-34 correspond to final results; note that col.26-34 correspond to size distribution (min, max, average) of different types of scaffolds (contaminant / containing at least 1 potential HGT / uncertain)

- 'settings': file with the date of the analyses & the versions of the softwares.
