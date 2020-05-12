import pandas as pd
import os 

##### load config and sample list #####

configfile: "config/config.yaml"

chromosomes = pd.read_table(config["chromosomes"], header=None)[0]

samples = pd.read_table(config["samples"])["sample"]

fastq1 = config["fastq1"] #path and format of fastq input 1
fastq2 = config["fastq2"] #path and format of fastq input 2

ref = config["ref"] #reference genome fasta (with index files in same directory)

outdir = config["outdir"] #base directory for output

##### target rule #####
rule all:
	input: 
		expand(os.path.join(outdir, "bams/{sample}.bam"), sample = samples), #Alignment output
		expand(os.path.join(outdir, "bams/{sample}_qc_report.html"), sample = samples) #QC output

# Align reads and filter
rule align:
	input: 
		ref,
		fastq1,
		fastq2
	output:
		os.path.join(outdir, "bams/{sample}.bam")
	conda:
		"envs/mapping.yaml"
	shell:
		"bwa mem -SP5 {input} | samblaster | samtools view -S -h -b -F 2316 > {output}"

#hard coded align:
#rule align:
#	input:
#		/u/home/j/jzou1115/project-zarlab/mm10/mm10.fa, #reference genome
#		/u/home/j/jzou1115/project-zarlab/pbc9/Hi-C_data/{sample}_L001_R1_001.fastq, #fastq file 1
#		/u/home/j/jzou1115/project-zarlab/pbc9/Hi-C_data/{sample}_L001_R1_002.fastq #fastq file 2
#	output:
#		/u/home/j/jzou1115/project-zarlab/HiC_pipeline_out/prelim/bams/{sample}.bam		
#	conda:
#		"envs/mapping.yaml"
#	shell:
#		"bwa mem -SP5 {input} | samblaster | samtools view -S -h -b -F 2316 > {output}"
		
#QC libraries
rule align_qc:
	input:
		os.path.join(outdir, "bams/{sample}.bam")
	params:
		prefix=os.path.join(outdir, "bams/{sample}")
	output:
		os.path.join(outdir, "bams/{sample}_qc_report.html"),
		report(os.path.join(outdir, "bams/{sample}_long.png"), caption = "report/align_qc.rst", category="QC Metrics", subcategory= "{sample}")
	conda:
		"envs/hic_qc.yaml"
	shell:
		"python scripts/hic_qc/hic_qc.py -b {input} -r -o {params.prefix}"


