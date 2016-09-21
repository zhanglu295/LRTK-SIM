SIMLR: A simulator for linked reads. The main programe is SIMLR.py, which is
wrote by python3
Prerequistion programs:
Compulsory:
1. dwgsim:short reads simulation 
Availability: http://github.com/nh13/DWGSIM
2. matplotlib: draw histogram
Availibility: http://matplotlib.org/
Option:
3. vcf-consensus: insert know variants to template sequence
Availability:https://vcftools.github.io/perl_module.html
4. tabix: index vcf files
Availability:https://sourceforge.net/projects/samtools/files/tabix/

All the installation path must be exported to the system path. Type dwgsim,
vcf-consensus, tabix to confirm no errors exist.

Command to run:
python3 SIMLR.py ./configure

The configuration fold includs parameter file for each library.
For example: two libraries are simulated

./configure
----------config1.txt (parameters for lib1)
----------config2.txt (parameters for lib2)

The simulated fastq files are written to the folder generated for each
library.
Two histograms are showed
1. Len_Molecule_hist.png: the length of molecules
number_Molecule_hist.png: the number of molecules in droplets
generated for



