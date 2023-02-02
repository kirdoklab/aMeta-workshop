---
title: "Fast alignment with Bowtie2"

---

## Bowtie2

The pipeline has a side branch for rapid analysis of the results based on an alignment with Bowtie2 and a post-mortem damage estimation with MapDamage. 
For this project, we won't be able to use the real Bowtie2 databases, because there were too heavy for the resources allocated for this workshop. So we are only going to look at the code but not run it and the output will be provided for you in the folder containing the expected results. 

## Bowtie2-build-l to build the index files

In order to run a Bowtie2 alignment, one needs a complete Bowtie2 database, in other words a .fna (fasta) file that has been indexed using the command bowtie2-build-l. This is the first part of the pipeline for the alignment step. You can therefore provide your own merged fna file for Bowtie2 to index if you wish to. However, if you are using one of our own custom databases, this step will be automatically ignored because its output already exists.

You can find the Bowtie2_build rule generating the index files below. It takes as input the path to the library.fna file that you have provided for the variable bowtie2_patho_db in your config.yaml file. For more information, refer to the main aMeta GitHub in the section about the config file. The rule then generates the index files of the Bowtie2 database **if they don't already exist**.

```
rule Bowtie2_Index:
    output:
        expand(
            f"{config['bowtie2_patho_db']}{{ext}}",
            ext=[
                ".1.bt2l",
                ".2.bt2l",
                ".3.bt2l",
                ".4.bt2l",
                ".rev.1.bt2l",
                ".rev.2.bt2l",
            ],
        ),
    input:
        ref=ancient(config["bowtie2_patho_db"]),
    log:
        f"{config['bowtie2_patho_db']}_BOWTIE2_BUILD.log",
    shell:
        ""bowtie2-build-l --threads {threads} {input.ref} {input.ref} > {log} 2>&1
```

Here is a simplified version of this code:

```
bowtie2-build-l --threads 1 resources/library.fna resources/library.fna > logs/Bowtie2_Build.log 2>&1
```

WARNING: No need to execute that line of code as we haven't downloaded the real databases to the server for this workshop and it would overwrite the index files! The trick for this rule is that the input file is declared as ancient with the "ancient" function of snakemake, making snakemake not rerun the rule if the output files already exists. But if you only run the bash command provided, you will overwrite the files.

### Regarding real life projects with the actual databases

For real life projects, these are the Bowtie2 databases we have made available for download:

+ Pathogenome: Bowtie2 index and helping files for following up on microbial pathogens
+ Bowtie2_Full_NT: Bowtie2 index for full NCBI NT (for quick follow up of prokaryotes and eukaryotes; also contains helping files for building the Malt database)

For more information and links to download the databases, please refer to the official GitHub of aMeta.

WARNING: if you are using the Bowtie2_Full_NT database, make sure that you have **unzipped** the files before running the pipeline, otherwise it will automatically overwrite the index files to regenerate them and will need a huge amount of memory and time to succeed, which you probably don't want to happen ;-)

## Bowtie2 alignment

```
rule Bowtie2_Pathogenome_Alignment:
    output:
        bam="results/BOWTIE2/{sample}/AlignedToPathogenome.bam",
        bai="results/BOWTIE2/{sample}/AlignedToPathogenome.bam.bai",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        db=rules.Bowtie2_Index.output,
    params:
        PATHO_DB=lambda wildcards, input: config["bowtie2_patho_db"],
    log:
        "logs/BOWTIE2/{sample}.log",
    shell:
        """bowtie2 --large-index -x {params.PATHO_DB} --end-to-end --threads 10 --very-sensitive -U {input.fastq} 2> {log} | samtools view -bS -q 1 -h -@ 10 - | samtools sort -@ 10 -o {output.bam} >> {log};"""
        """samtools index {output.bam}"""
```

Here is a simplified version of this code:

```
bowtie2 --large-index -x resources/library.fna --end-to-end --threads 10 --very-sensitive -U results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz 2> logs/BOWTIE2/{sample}.log | samtools view -bS -q 1 -h -@ 10 - | samtools sort - -@ 10 -o results/BOWTIE2/{sample}/AlignedToPathogenome.bam >> logs/BOWTIE2/{sample}.log
samtools index results/BOWTIE2/{sample}/AlignedToPathogenome.bam
```

Again, no need to run that line of code. You can find how the output is supposed to look like in the folder with expected results.
Basically, this line asks Bowtie2 to align each fastq file to the Bowtie2 database. This will generate a sam file containing the information of alignment for each DNA sequence and the reference genome to which it aligns to. The sam file is then directly changed into a bam file using samtools, sorted and later indexed so that it is ready to be used. In order to run this line, you would still need to replace {sample} with the name of each sample that you have provided in the samples.tsv file and run this line for each of them.
