README FOR DELSCAN 0.1
AUTHOR: Don Conrad 
DATE: 20/10/2010


########### INTRO ###############

What this program does: Given a set of SNP genotypes from a parent-offspring trio, delscan uses a dynamic programming algorithm to identify stretches of genome where the offspring and one parent are homozygous for all SNP genotype calls, which may be indicative of a transmitted deletion. To increase the specificity of this approach it further requires the presence of Mendelian errors (MEs) in order to declare a candidate region, and these MEs must all be consistent with the deletion event (either maternal or paternal transmission) suggested by the patterns of homozygosity. 

What this program does not do: There is no formal measure of confidence associated with the regions reported, and the results are meant for exploratory use only. In general, the number of Mendelian errors within a deletion interval is the strongest indicator of the probability that this deletion is real. By default the program will only report intervals with 2 or more MEs. In an analysis of HapMap data, completed in 2005 and published in 2006, I estimated the False Discovery Rate of the method, when requiring 2 or more MEs, to be 10% - 15%, depending on the population being analyzed. I suspect that due to improvements in array chemistry and genotyping algorithms, the FDR would be even lower on contemporary data. 


In this archive you will find the following:

Code:

run_dscan.pl -  a wrapper for delscan that takes a set of three "23andMe" formatted genotype files, creates single-chromosome input files for delscan, calls delscan and creates a single merged output file.

delscan.c - source code for the main program. 

Binaries:

delscan_darwin_64 - 64bit binary compiled with i686-apple-darwin8-g++-4.0.1
delscan_linux_64  - 64bit binary compiled with g++ 4.3.2, in Debian Linux environment (4.3.2, "Lenny"). 
delscan_linux_32  - 32bit binary compiled with g++ in Debian Linux environment ("Etch")


Example input files for a new project, using data from the International HapMap Consortium (www.hapmap.org):

hapmap.fam

indat/NA10863_HM3.txt
indat/NA12234_HM3.txt
indat/NA12264_HM3.txt

###### BEFORE YOU START #########

If you are using a Linux or Mac machine then you may be able to use the binaries provided above. Just copy 
or rename the binary of choice to "delscan". Otherwise you will need to compile the code from scratch

Compiling code
--------------

you'll need a C++ compiler such as g++.

from the command line:
g++ -o delscan delscan.c -lm



####### Using "run_dscan.pl"  #######

USAGE : perl run_dscan.pl <input directory> <family file>

This PERL wrapper opens up <family file> and reads in the names of the three genotype datasets (one for mom, dad, child), which should be located in the <input directory>. To minimize memory usage it splits the input data into temporary single chromosome datafiles. run_dscan.pl then calls the C program and creates a single output file. 


"Family File"
This should be a simple text file with three rows as in the following:

MOTHER name_of_mothers_data_file MOTHER_ID
FATHER name_of_fathers_data_file FATHER_ID
CHILD name_of_childs_data_file	CHILD_ID

The first column should always contain these labels ("MOTHER","FATHER","CHILD") in that order. The other columns are variable values that must be supplied by the user.

EXAMPLE: I've included some example input files with this package so that you can test that code is running correctly, and practice running the wrapper. The files are 

hapmap.fam (example family file)

A folder, "indat", containing genotype data on chromosome 8 from a HapMap trio:

NA10863_HM3.txt 
NA12234_HM3.txt
NA12264_HM3.txt

To test that everything is working, run 

run_dscan.pl indat hapmap.fam

This should produce a file called "hapmap.fam.dels" containing 6 deletion intervals

####### WRAPPER OUTPUT ##########

The make_input.pl creates a single output file. The root name of this output is the name of the family file, with the suffix ".dels". If my family file is named "conrad.fam" then the output is "conrad.fam.dels". There will be one row per potential deletion. If there is no output in the file, and everything ran properly, then no deletions were found.

The format is as follows:

chromosome start end nsnp length nme carrier offspring_id end_flag

Output Field Descriptions
-------------------------

start, end - the edges of the deletion are defined by the first SNP incompatible with the deletion on either side. The breakpoints 
reported are INCLUSIVE of these first incompatible SNPs, and thus represent the maximum possible extent of the deletion. It is useful 
to report breaks this way, as the inter-snp density on some SNP chips may be large with respect to the size of variants detected by the program.

nsnp         - number of SNPs within the deletion stretch reported, again this include 2 sites that are not within the deletion

length       - length of deletion stretch in basepairs

nme          - number of mendelian errors in the stretch

carrier      - which parent transmitted the variant to the child; 1 = mother 2= father.

offspring_id - sample id of trio offspring, not really useful unless running multiple trios


####### RUNNING DELSCAN WITHOUT WRAPPER ##########

Usage:

./delscan genotype_file pedigree_file output_format input_format chromosome_number

Pedigree file is a standard pedfile format. Each row corresponds to a family member. The fields for each row are

Family_ID Sample_ID Fathers_ID Mothers_ID Sex 

The value of the "Family_ID" doesn't really matter, as long as everyone from the same family has the same Family_ID. The Sample IDs must be identical to the IDs that appear in the genotype file.

Genotype file formats supported: 

0 HapMap format 
1 PedCNV
2 BEAGLE format

