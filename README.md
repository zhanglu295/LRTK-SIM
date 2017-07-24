# SIMLR: Linked reads simulator for 10X Chromium. 
## Prerequired software

[python3](https://www.python.org/downloads/source/) and [matplotlib](https://matplotlib.org/) are requeired.
##facilitated program
`gen_fasta.py` can generate two haploid fasta files with variants from vcf and reference genome sequence.

`usage: python gen_fasta.py -v sample.vcf -r ref.fasta -p newref -o ./work`

`-v --vcf, the input path of compressed or uncompressed vcf file`

`-r --reference, the path of compressed or uncompressed ref file`

`-p --prefix, prefix of new reference files`

`-o --out, the path to output`

`-h --help, help info`


## Basic usage

`python SIMLR.py <path to config files>`

Multiple libraries are allowed to simulated simutaneously and the parameters are included in the correspinding config files. The config files must be named as libX, where X should be a number as the library id.

E.g. `python SIMLR.py ./configure`

The configure folder includes two config files for two different libraries:

config1.txt (for lib1) and config2.txt (for lib2)

The simulated fastq files are written to the folder generated for each
library.

## Important parameters in config file

Two examples of config file (config1.txt and config2.txt) are prepared in configure folder.

`line2` and `line3`: `Path_Fastahap1` and `Path_Fastahap2`, the target genomic sequences to be simulated. At most, SIMLR allows two fasta files to represent two haplotypes. The diploid genomic sequences can be generated by [vcf-consensus](https://vcftools.github.io/perl_module.html) to insert variants to the reference genome. If only would like to perform haploid genome simulation, just remove `Path_Fastahap2=XXX` in `line3` 

`line5`: `processors`, the maximum number of CPUs are allowed to use

`line7`: `CF`, coverage of long fragments

`line9`: `CR`, covergae of short reads

`line11`: `N_FP`, the average number of molecules for each droplet

`line13`: `Mu_F`, the average length for long fragment (Kb)

`line15`: `SR`, length of short reads (bp)

`line21`: `Error_rate`, sequencing error rate

`line27`: `Mu_IS`, the average of insert size for short reads (bp)

`line29`: `Std_IS`, standard deviation of insert size for short reads (bp)

`line33`: `Hap`, Haploid (Hap=1) or Diploid (Hap=2)
## Output
`SIMLR` generates one folder for each library named after libX (X is the library id) in the path of config files. Two compressed fastq files, started by read-I1 (sample index) and read-RA (short reads) are generated, which are required by 10X official pipelines [Long Ranger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) and [Supernova](https://support.10xgenomics.com/de-novo-assembly/software/overview/welcome).

Also two histograms are generated

1. Len_Molecule_hist.png: the Molecule length distribution 

2. Number_Molecule_hist.png: the distribution of molecule number in droplets

Contact: zhanglu295@gmail.com
