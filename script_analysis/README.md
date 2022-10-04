# Example real data analysis pipeline:

In the CAPG manuscript, we collected WGS data from fourteen peanut accessions, aligned them to the [Tifrunner reference](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_003086295.2/), and genotyped 1,000 selected target regions.
Here, we demonstrate a small portion of the analysis. For this example, chr 1 and chr 11 of the [Tifrunner assembly](https://peanutbase.org/data/v2/Arachis/hypogaea/genomes/Tifrunner.gnm2.J5K5/) is used as refA.fa and refB.fa. 
 
All commands are run from the top directory of the github repository.
**Unfortunately, the following demonstration does not work because of a disagreement between the Tifrunner assembly and the version we used. We are working on it.**

- Big data files are not stored on the github repository, so the first step is to download the files you will need to run this example.

	- Download the fastq file using `wget`
	```
	wget https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR4124062&display=download
	```

        OR the SRA toolkit
        ```
        fasterq-dump --split-files SRR4124062
        ```
- Install the following packages:
	- bwa2 or an aligner of your choice
	- samtools
	- bioawk
	- Mummer4

- Perform alignment of example fastq with reference A and B genome using `bwa2 mem` (or using your choice of aligner)
```
bwa2 index data/peanut/refA.fa
bwa2 index data/peanut/refB.fa
bwa2 mem data/peanut/refA.fa data/peanut/SRR4124062_1.fastq data/peanut/SRR4124062_2.fastq | samtools sort -o data/peanut/sam/SRR4124062_A.bam
bwa2 mem data/peanut/refB.fa data/peanut/SRR4124062_1.fastq data/peanut/SRR4124062_2.fastq | samtools sort -o data/peanut/sam/SRR4124062_B.bam 
```

- Generate reference sam file
```
nucmer --sam-long=peanut_ref --mum data/peanut/refA.fa data/peanut/refB.fa
bioawk -c fastx '{ print "@SQ " "SN:" $name, "LN:" length($seq) }' < data/peanut/refA.fa > header.txt
cat data/peanut/header.txt data/peanut/peanut_ref.sam
```

- Subset the sam files according to the target region of interest (10 target regions have been selected in targets.txt for this demonstartion)
```
cat /data/peanut/targets.txt | awk -F "[\t : -]" '{ printf "%s\t%s\t%s\n", $2,$3,$4}' > /data/peanut/subset_A.bed
cat /data/peanut/targets.txt | awk -F "[\t : -]" '{ printf "%s\t%s\t%s\n", $5,$6,$7}' > /data/peanut/subset_B.bed
samtools view -h -L /data/peanut/subset_A.bed > data/peanut/sam/SRR4124062_A_subset.bam
samtools view -h -L /data/peanut/subset_B.bed > data/peanut/sam/SRR4124062_B_subset.bam
```

- Now you are ready to run the genotyping pipeline.
There is nothing pretty about this pipeline, but it gets the job done for now.
```
script_analysis/genotype_WGS.sh		# run capg_wgs: vcf files in data/peanut/vcf
script_analysis/err2info.sh		# convert captured stderr output of capg_wgs to "info" format
script_analysis/combine_infos.sh	# combine "info" files
script_analysis/get_metrics.R		# compute metrics from "info" files
```
- The metrics for heterozygote genotyping are in the file called `results/peanut/CAPG_PL_peanut_het.txt`.
The metrics for allelic (homologous) and homoeologous SNP calling are in the file called `results/peanut/CAPG_PL_peanut_ho.txt`.

# Example simulation analysis pipeline:

In the CAPG manuscript, we simulated data from 50 individuals over seven conditions.
Here, we demonstrate a simulation with 5 individuals over 12 conditions.

- CAPG pipeline: You must have succeeded in compiling `capg_sim`, which requires the R math library and installed `bwa` and [art_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
```
./script_analysis/simulate_data.sh		# calls capg_sim
./script_analysis/capg_genotype_simulation.sh	# calls capg_wgs
./script_analysis/err2info_simulation.sh	# calls err2info_simulation.py
./script_analysis/get_metrics_simulation.R	# computes metrics
```
VCF files are placed in `data/simulation/homrHRmmMM/covCVG/vcf`, where `HR` is homoeologous rate (0.005 or 0.007), `MM` is subgenome reference mismatch rate (0.000, 0.001, 0.010), CVG is read coverage (5 or 10 per chromosome).
Computed metrics for every site are in `data/simulation/results/CAPG_PL_HR_CVG_mmMM_het.txt` for heterozygous genotype calls or `data/simulation/results/CAPG_PL_HR_CVG_mmMM_het.txt` for SNP calls (allelic or homoeologous).
- GATK pipeline: You must have simulated the data using the above CAPG pipeline and installed `samtools`, `bwa`, and `gatk`.
```
./script_analysis/simulation_GATK.sh		# creates vcf files
./script_analysis/vcf2info_simulation.sh	# calls vcf2info.py
./script_analysis/get_metrics_GATK_simulation.R	# computes metrics
```
VCF files are placed in `data/simulation/homrHRmmMM/covCVG/GATK/vcf` for given `HR`, `MM`, and `CVG` value (see above).
Computed metrics are in `data/simulation/results/GATK_PL_HR_CVG_mmMM_het.txt` or `data/simulation/results/GATK_PL_HR_CVG_mmMM_hom.txt`.
