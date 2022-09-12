# CAPG - Comprehensive Allopolyploid Genotyper

This software genotypes targeted genomic regions in allotetraploids using whole genome sequencing (WGS) reads aligned to both reference subgenomes.
(For a version that handles targeted amplicons, see [capg_amp](#amplicon).)
The main genotyper is written in C as a standalone executable.

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Installation](#installation)
1. [Required Inputs](#input)
1. [Output](#output)
1. [Other Command-Line Options](#options)
1. [Tutorial](#tutorial)
1. [How to Cite](#cite)
1. [Contact](#contact)
1. [Amplicon Version](#amplicon)

# Prerequisites <a name = "prerequisites" />

The genotyper ```capg_wgs``` requires the [C compiler from the GNU Compiler Collection (GCC)](https://gcc.gnu.org/), [CMake](https://cmake.org/), and [Samtools](https://www.htslib.org/download/) executable installed on your system.

There is also a data simulator ```capg_sim``` available.
It additionally requires the [R Standalone Math Library](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library).  Often, the Rmath library (```libRmath.a``` or ```libRmath.so``` for Linux or ```libRmath.dylib``` for MacOS) will be installed with R, but not always.  Here are some other locations for the library.

- r-mathlib on [Ubuntu](https://ubuntu.com/) and [Debian](https://www.debian.org/)
- libRmath on [Fedora](https://ubuntu.com/), [CentOS](https://centos.org/), [Mageia](https://www.mageia.org/en/), and [Mandriva](https://www.openmandriva.org/)
- Or if all else fails, you can install the Rmath standalone library from the repository [https://github.com/statslabs/rmath](https://github.com/statslabs/rmath)

If RMathLib is not installed on your system, everything should be fine except ```capg_sim``` will not be compiled.

# Installation <a name = "installation" />

1. Clone the repository.

    ```sh
    git clone https://github.com/Kkulkarni1/CAPG.git
    ```

2. Compile CAPG. The executable is called ```capg_wgs```.  It will appear in the ```CAPG/src/wgs``` directory.

   ```sh
   cd CAPG/src/wgs
   cmake .
   make
   ```

3. Install CAPG. Copy the executable to wherever you need it. If you have root privileges, you can install it into the system path, for example:

   ```sh
   sudo cp capg_wgs /usr/local/bin
   ```
   
# Required Inputs <a name="input" />

The software requires multiple input files.

1. Fasta files containing subgenomic references, one subgenome per file. Pass them in via the ```--fsa_files``` command-line option.

2. SAM files containing the reads aligned to each subgenomic reference separately.  Pass them in via the ```--sam_files``` command-line option.

3. SAM file containing the alignments of selected target regions in each subgenome to each other.  Pass it in via the ```--geno``` command-line option.

It also requires one command-line option.

You must name the target regions to genotype, including the chromosome name and the start and end positions relative to the whole chromosome.
For example, ```chr1:1-10``` means you want to genotype from position 1 to 10 (1 based) in chromosome 1. Use ':' to seperate chromosome name and region index. Use '-' to seperate start and end positions.
Pass these in by the ```--ref_names``` command-line option.

We have used [MUMmer4](https://github.com/mummer4/mummer) to produce the SAM file for the alignment of the targeted region(s).
For example, the following command will output a SAM file called ```ref.sam```.

```
nucmer --sam-long=ref --mum target_A.fa target_B.fa
```


# Output <a name="output" />

The genotyping output for each subgenome are stored in [VCF files](https://samtools.github.io/hts-specs/VCFv4.2.pdf), one per subgenome, if the `--vcf_files` command-line option is used.
The name used to identify the current individual in the VCF file output can be provided with the `--name` option.
In addition, the program will extract the target regions into FASTA files, by default called `extracted0.fsa` and `extracted1.fsa`, though you can change the prefix `extracted` with the `-j` command-line option.
These files are currently deleted unless the program terminates unexpectedly, so you can place these in a temporary directory.
The command also currently produces //a lot// of output to `stderr` that you may wish to capture and examine.


# Command-Line Options <a name = "options" />

Please run `./capg_wgs -h` for detailed information about all available options.


# Tutorial <a name = "tutorial" />

All the files used and created in this tutorial are in the `data` folder. 
In this example we will genotype positions 1 to 5000 of both subgenomes assuming they are homoeologous.
Finally, we store the output in VCF files, whose names are provided via the `--vcf_files` option.

From the `src/wgs` directory, the command line for genotyping is:

```
./capg_wgs --ref_names Genome_A:1-5000 Genome_B:1-5000 --sam_files ../../data/aln0A.sam ../../data/aln0B.sam --fsa_files ../../data/refA.fa ../../data/refB.fa --geno ../../data/ref.sam -equal --vcf_files ../../data/A.vcf ../../data/B.vcf
```

The positions with no coverage in the first genome will not be outputed.

More detailed tutorials demonstrating real data analysis of peanut data and an extensive simulation, can be found [here](https://github.com/Kkulkarni1/CAPG/tree/main/script_analysis).

# How to Cite <a name = "cite" />

- This work is under review.  Please see [bioarxiv](https://biorxiv.org/cgi/content/short/2022.04.21.488948v1).

# Contact <a name = "contact" />

If you have any problems with this software, please contact:

Roshan Kulkarni (roshank@iastate.edu) or Karin S. Dorman (kdorman@iastate.edu)

# Amplicon Version <a name = "amplicon" />

We also have a similar software for genotyping amplicon sequences, here we briefly mention how to use it.

## Installation

Compile CAPG for amplicon. The executable is called ```capg_amp```.  It will appear in the ```CAPG/src/amplicon``` directory.

   ```sh
   cd CAPG/src/amplicon
   make capg_amp
   ```
## Command-Line Options

```
CAPG_AMP(1)

NAME
	capg_amp - genotype tetraploids

SYNOPSIS
	capg_amp --sam_files SAM1 SAM2 --fasta_files FSA1 FSA2 --ref_names REF1 REF2
		[[--genotype_by_clustering [--alignment FILE1 FILE2]]
		[--sample INT --min-subgenomic-coverage FLOAT]
		[--min INT --max INT --expected-errors FLOAT --indel INT --loglik FLOAT
		 --min-posterior FLOAT --secondary --soft-clipped INT]
		[--coverage FLOAT --biallelic FLOAT --equal_coverage_test [FLOAT1 FLOAT2]]
		[--drop INT --amplici EXE [--amplici-f FILE --amplici-o STRING --amplici-l FLOAT]]
		[--error_file|--error_data FILE] ...]

DESCRIPTION
	capg_amp genotypes allotetraploids using reads in SAM1 and SAM2 aligned to
	REF1 and REF2 references from fasta files FSA1 FSA2.
	SAM files typically contain reads from a single individual, genotype, or
	accession aligned to multiple amplified targets, but capg_amp genotypes
	one individual at one amplicon.

OPTIONS

Input (required):
	--fasta_files FILE1 FILE2
		Subgenomic reference fasta files (Default: none).
		DEPRECATED: see --ref_fasta_files
	--ref_fasta_files FILE1 FILE2
		Subgenomic reference fasta files (Default: none).
	--sam_files FILE1 FILE2
		SAM files with reads aligned to each subgenome (Default: none)
	--ref_names STRING1 STRING2
		Names of subgenomic reference target regions (Default: none)
	--ref_alignment FILE
		SAM file containing alignment of references (Default: none)

Output (optional):
	--display_alignment
		Display alignments in stderr output (Default: no).
	--vcf_files FILE1 FILE2
		Genotyping output in one vcf file per subgenome (Default: none).
	--subref_fasta_files FILE|FILE1 FILE2
		Subsetted reference regions output to these FASTA files (Default: subsetted_refs[12].fsa)
	--gl
		Toggle GL output to vcf files (Default: yes).
	--name STRING
		Name of accession/individual/genotype; used in vcf header (Default: sample).

Estimation/Inference (optional):
	--genotype_by_clustering
		Genotype by clustering (Default: no).
		Requires command-line arguments --amplici and --clustalo or --mafft.
	--clustalo FILE
		The clustal omega executable (Default: (null)).
	--mafft FILE
		The mafft executable (Default: (null)).
	--alignment FILE1 FILE2
		Alignment input FILE1 and output FILE2 (Default: selected_haplotypes.fa selected_haplotypes.co.fa)
	--misalignment_rate FLOAT
		Maximum allowed subgenomic misalignment rate in [0, 0.5) (Default: 0.30).
		Tolerated proportion of reads from one subgenome aligning to the other.

Screening paralogs and other contaminants (optional):
	-p, --drop INT
		Drop reads aligning to paralogs; -1 to automate (Default: -1)
		Drop INT of four most abundant haplotypes if specified.
		Requires command-line argument --amplici.
	--amplici EXE
		The amplicon denoiser software (Default: none).
		Writes auxiliary files "amplici.fastq" "amplici.fa", and "amplici.out".
		See https://github.com/DormanLab/AmpliCI for more information.
	--amplici_fastq FILE
		Selected reads output to this FASTQ file for denoising (Default: amplici.fastq).
	--write-fastq [FILE]
		Write fastq file for AmpliCI and quit (Default: no).
		See --amplici_fastq to name the file.
		With optional argument write named fastq file after AmpliCI paralog filtering (Default: none).
```
	
