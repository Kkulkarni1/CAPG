# CAPG - Comprehensive Allopolyploid Genotyper

This software is implemented for genotyping alleotetraploids' targeted genome regions using screened reads aligned to the whole genome references.
It has primarily been designed as a standalone executable.
Skip to [Installation](#installation) to see how.

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Installation](#installation)
1. [Command-Line Options](#options)
1. [Input Files](#input)
1. [Output](#output)
1. [Troubleshooting](#troubleshooting)
1. [Tutorial](#tutorial)
1. [How to Cite](#cite)
1. [Acknowledgements](#acknowledgements)
1. [Contact](#contact)

# Prerequisites <a name = "prerequisites" />
- Rmath, the [R Standalone Math Library](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library).  Often, the Rmath library (libRmath.a or libRmath.so for Linux or libRmath.dylib for MacOS) will be installed with R, but not always.  Here are some other locations for the library.
	- r-mathlib on [Ubuntu](https://ubuntu.com/) and [Debian](https://www.debian.org/)
	- libRmath on [Fedora](https://ubuntu.com/), [CentOS](https://centos.org/), [Mageia](https://www.mageia.org/en/), and [Mandriva](https://www.openmandriva.org/)
	- Or if all else fails, you can install the Rmath standalone library from the repository [https://github.com/statslabs/rmath](https://github.com/statslabs/rmath)

- [Samtools](http://www.htslib.org/download/) should be installed into the sysmtem path /usr/local/bin.

# Installation <a name = "installation" />

1. Clone the repository.

    ```sh
    git clone https://github.com/Kkulkarni1/CAPG.git
    ```

2. Compile CAPG. The executable is called ```capg_wgs```.  It will appear in the ```CAPG/ksd-roshan_wgs``` directory you are currently in.

   ```sh
   cd CAPG/ksd-roshan_wgs
   make capg_wgs
   ```

3. Install CAPG. Copy the executable to wherever you need it. If you have root privileges, you can install it into the system path, for example:

   ```sh
   sudo cp capg_wgs /usr/local/bin
   ```
   
# Command-Line Options <a name = "options" />

Please run `./capg_wgs -h` for detailed information about all available options.  That help is repeated below.
```
CAPG_WGS(1)

NAME
	capg_wgs - genotype tetraploids

SYNOPSIS
	capg_wgs --sam_files <fsam1> <fsam2> --fsa_files <fsa1> <fsa2> --ref_names <sref1> <sref2> --g <refsam> --j <reffsa>
		[--sample <int> --min-subgenomic-coverage <dbl>]
		[--min <int> --max <int> --expected-errors <dbl> --indel <int> --loglike <dbl> --secondary --soft-clipped <int> --coverage <dbl>]
		[ --o <fout> --error_file|--error_data <ferr>] ...]

DESCRIPTION
	capg_wgs genotypes tetraploids' targeted genome regions using screened reads in <fsam1> and <fsam2> aligned to the whole genome references from fasta files <fsa1> <fsa2>

NOTICE
	capg_wgs requires the reference names (appeared in <refsam>) to contain start (0 based) and end position (1 based) of the targeted genome regions relative to the whole genome. Using ':' to seperate original genome name (appeared in <fsam1> and <fsam2>) and region index, '-' to seperate start and end position (e.g. chr1:0-11, which starts at 1st position in 'chr1' and end at 11th position, length is 10 bases)

OPTIONS
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Input: all required
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	--fsa_files <fsa1> <fsa2>
		Specify fasta files containing subgenomic reference sequences [use samtools faidx to index the files] (Default: none).
	--sam_files <fsam1> <fsam2>
		Specify sam files containing alignments (Default: none)
	--ref_names <sref1> <sref2>
		Specify names of subgenomic references for target region; must exist in sam  files (Default: none)
	--geno <refsam>
		Specify name of sam file of aligning <sref2> to <sref1> (Default: none)
	--j <reffsa>
		Specify name of targeted fsa files to be extracted  by samtools [Set samtools in system PATH] (Default: none)
	+++++++++++++++++
	Output:  optional
	+++++++++++++++++
	--display_alignment
		Display alignments in stderr output (Default: no).
	--vcf_files FILE1 FILE2
		Genotyping output in one vcf file per subgenome (Default: none).
	--gl
		Toggle GL output to vcf files (Default: no).
	--name STRING
		Name of accession/individual/genotype; used in vcf header (Default: none).
	--o <fout>
		Specify file with the final output (Default: none).
	++++++++++++++++++++++++++++++
	Error recalibration:  optional
	++++++++++++++++++++++++++++++
	--error_file <ferr>
		Specify file with estimates of error rates (Default: none).
	--error_data <ferr>
		Specify file to output observed errors (Default: none).
	+++++++++++++++++++++++++++++++++++++++++++
	Screening reads, coverage checks:  optional
	+++++++++++++++++++++++++++++++++++++++++++
	--po <pint>
		Equal coverage test. (Default: no)
	--expected_errors <dbl>
		Discard reads with more than <dbl> expected errors (Default: inf).
	--indel <i>
		Drop reads with alignments containing more than <i> indels (Default: 2147483647)
	--loglik <l>
		Drop reads with log likelihood less than <l> (Default: -inf)
	--biallelic FLOAT
		Skip site if third allele >100*FLOAT% of minimum subgenomic coverage (Default: 0.5).
	--min <dbl>
		Drop reads shorter than <dbl> (Default: 0)
	--max <dbl>
		Drop reads longer than <dbl> (Default: 2147483647)
	--secondary
		Drop secondary alignments (Default: yes)
	--coverage <c>
		Tuning parameter for penalty (Default: 1.0).
	--eq
		Post-hoc test of equal coverage of homologous chromosomes. (Default: yes)
	--soft-clipped <s>
		Drop reads where either alignment is clipped by <s> or more nucleotides (Default: 2147483647)
	--unmapped
		Drop reads unmapped in either alignment (Default: yes)
	+++++++++++++++++++++++++++++++
	Estimation/Inference:  optional
	+++++++++++++++++++++++++++++++
	--min_subgenomic_coverage <c>
		Abort if subgenomic coverage drops below <c> (Default: 5.0).
```

# Input Files <a name="input" />

The software requires multiple input files.

1. Fasta files containing subgenomic references.

2. SAM files containing the reads alignments to each of the references.

3. SAM file containing the alignments of selected regions for genotyping. Our software has a requirement of the names of the genotyping regions, they needs to contain start and end position relative to the whole genome and the genome name. Using ':' to seperate original genome name and region index, '-' to seperate start and end position. For example, chr1:0-10 means we will genotype from position 1 to 10 (1 based).

We recommand using [MUMmer4](https://github.com/mummer4/mummer) to produce the SAM file for selected regions, command is given below.

	```
	nucmer --sam-long=combine --mum target_A_homeolog.fsa target_B_homeolog.fsa
	```
