# Example real data analysis pipeline:

In the CAPG manuscript, we collected WGS data from fourteen peanut accessions, aligned them to the [Tifrunner reference](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_003086295.2/), and genotyped 1,000 selected target regions.
Here, we demonstrate a small portion of the analysis.
All commands are run from the top directory of the github repository.
**Unfortunately, the following demonstration does not work because of a disagreement between the Tifrunner assembly and the version we used. We are working on it.**

- Big data files are not stored on the github repository, so the first step is to download the files you will need to run this example.

	- Download the [Tifrunner assembly](https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/GCF_003086295.2/download?filename=GCF_003086295.2.zip) and overwrite the github stubs for the subgenomic reference sequences:
	```
	wget https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/GCF_003086295.1/download?filename=GCF_003086295.1.zip GCF_003086295.1.zip
	cat ncbi_dataset/data/GCF_003086295.2/chrArahy.0[1-9].fna ncbi_dataset/data/GCF_003086295.2/chrArahy.10.fna | awk -f script_analysis/chg_names.awk > data/peanut/tet_A.fa	# overwrite github stub
	cat ncbi_dataset/data/GCF_003086295.2/chrArahy.1[1-9].fna ncbi_dataset/data/GCF_003086295.2/chrArahy.20.fna | awk -f script_analysis/chg_names.awk > data/peanut/tet_B.fa	# overwrite github stub
	```
	- Download the target-aligned SRR4124062 reads from [OSF storage](https://osf.io/uezgp/files/osfstorage):
	```
	wget https://osf.io/download/631d4476db9397378e11f644/ data/peanut/sam/SRR4124062_A.subset.bam
	wget https://osf.io/download/631d44929a7d513523903ae0/ data/peanut/sam/SRR4124062_B.subset.bam
	samtools view -h data/peanut/sam/SRR4124062_A.subset.bam data/peanut/sam/SRR4124062_A.subset.sam	# overwrite github stub
	samtools view -h data/peanut/sam/SRR4124062_B.subset.bam data/peanut/sam/SRR4124062_B.subset.sam	# overwrite github stub
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

- CAPG pipeline
```
./script_analysis/simulate_data.sh		# calls capg_sim
./script_analysis/capg_genotype_simulation.sh	# calls capg_wgs
./script_analysis/err2info_simulation.sh	# calls err2info_simulation.py
./script_analysis/get_metrics_simulation.R	# computes metrics
```
VCF files are placed in `data/simulation/homrHRmmMM/covCVG/vcf`, where `HR` is homoeologous rate (0.005 or 0.007), `MM` is subgenome reference mismatch rate (0.000, 0.001, 0.010), CVG is read coverage (5 or 10 per chromosome).
Computed metrics for every site are in `data/simulation/results/CAPG_PL_HR_CVG_mmMM_het.txt` for heterozygous genotype calls or `data/simulation/results/CAPG_PL_HR_CVG_mmMM_het.txt` for SNP calls (allelic or homoeologous).
- GATK pipeline: you must have installed `samtools`, `bwa`, and `gatk`.
```
./script_analysis/simulation_GATK.sh		# creates vcf files
./script_analysis/vcf2info_simulation.sh	# calls vcf2info.py
./script_analysis/get_metrics_GATK_simulation.R	# computes metrics
```
VCF files are placed in `data/simulation/homrHRmmMM/covCVG/GATK/vcf` for given `HR`, `MM`, and `CVG` value (see above).
Computed metrics are in `data/simulation/results/GATK_PL_HR_CVG_mmMM_het.txt` or `data/simulation/results/GATK_PL_HR_CVG_mmMM_hom.txt`.
