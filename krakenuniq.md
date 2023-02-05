---
title: "Prescreening with KrakenUniq"
---

## KrakenUniq prescreening step

The first part of the pipeline includes a prescreening step with KrakenUniq. It is not an aligner but a k-mer based classifier. It means that it has a database with information about unique regions of the DNA of each organism and classifies the reads by comparing them to this database. It can assign the reads to different taxonomic clade depending on their uniqueness. If a read is specific to a species, it will be assigned to that species, if it is specific to a genus, it will be assigned to that genus and so on. 

KrakenUniq is fast and makes it possible to screen aDNA samples against a **database as large as possible**. With an aligner, this is rarely possible because it would take too much compute power and time to align to a very large database. The KrakenUniq paper suggests that KrakenUniq is not less accurate compared to alignment tools such as `BLAST` and `MEGAN`. However, we still want to have a secondary verification with an aligner later on in the pipeline, since we need to align the reads in order to verify the quality of the alignment and generate statistics. 

## KrakenUniq database

We have made two KrakenUniq databases available through the SciLifeLab repository. Description of the databases and information on how to download them can be found on the [official aMeta GitHub](https://github.com/NBISweden/aMeta).

## Classification with KrakenUniq

**Please note that it will not be possible to run this first step during this workshop due to issues with storage space.** If you followed the course setup correctly, you should have the expected output files of this rule in your `results` folder. 

```
rule KrakenUniq:
    output:
        report="results/KRAKENUNIQ/{sample}/krakenuniq.output",
        seqs="results/KRAKENUNIQ/{sample}/sequences.krakenuniq",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
    params:
        DB=config["krakenuniq_db"],
    log:
        "logs/KRAKENUNIQ/{sample}.log",
    shell:
        "krakenuniq --preload --db {params.DB} --fastq-input {input.fastq} --threads {threads} --output {output.seqs} --report-file {output.report} --gzip-compressed --only-classified-out &> {log}"
```

Here is a shell version of this code:

```bash
for sample in $(ls results/CUTADAPT_ADAPTER_TRIMMING/*.fastq.gz); do
       sample_name=$(basename $sample .fastq.gz)
       krakenuniq --preload --db $DBNAME --fastq-input ${sample_name} --threads 4 --output ${sample_name}.sequences.krakenuniq --report-file ${sample_name}.krakenuniq.output --gzip-compressed --only-classified-out &> logs/KRAKENUNIQ/${sample_name}.log
done
```
In summary, this command loops through the sample files created by Cutadapt and classifies the reads using KrakenUniq.

## Filter the KrakenUniq output

After running KrakenUniq, we need to filter its output to remove as many false positives as possible. In order to do so, we select the species taxonomic level for each organism and we filter the results according to the amount of kmers and the amount of taxReads (reads specific to a species). A suggested value for this parameters is 1000 unique kmers to make sure that at least 1000 unique regions of the organim are covered and 200 taxReads in order to have enough reads to verify if the organism is ancient after alignment. 

```
rule Filter_KrakenUniq_Output:
    output:
        filtered="results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered",
        pathogens="results/KRAKENUNIQ/{sample}/krakenuniq.output.pathogens",
        pathogen_tax_id="results/KRAKENUNIQ/{sample}`taxID.pathogens",
    input:
        krakenuniq="results/KRAKENUNIQ/{sample}/krakenuniq.output",
        pathogenomesFound=config["pathogenomesFound"],
    log:
        "logs/FILTER_KRAKENUNIQ_OUTPUT/{sample}.log",
    params:
        exe=WORKFLOW_DIR / "scripts/filter_krakenuniq.py",
        n_unique_kmers=config["n_unique_kmers"],
        n_tax_reads=config["n_tax_reads"],
    shell:
        """{params.exe} {input.krakenuniq} {params.n_unique_kmers} {params.n_tax_reads} {input.pathogenomesFound} &> {log}; """
        """cut -f7 {output.pathogens} | tail -n +2 > {output.pathogen_tax_id}"""
        
```

Here is a shell version of this code:

```bash
# Suggested parameters for the number of unique kmers and the number of taxReads
n_unique_kmers=1000
n_tax_reads=200

# Create the log folder
mkdir -p logs/FILTER_KRAKENUNIQ_OUTPUT

# Loop through the samples and filter the KrakenUniq output according to three thresholds. It should have at least 1000 unique kmers and 200 reads and it should be at the Species level (not Genus, not Family, subspecies or else). This part is implemented in the python script.
# The second command extracts the column containing the taxID information for the species that match a pathogen in the pathogenFound.very_inclusive.tab.
for sample in $(ls results/CUTADAPT_ADAPTER_TRIMMING/*.fastq.gz); do
        sample_name=$(basename $sample .trimmed.fastq.gz)
        python scripts/filter_krakenuniq.py results/KRAKENUNIQ/${sample_name}/krakenuniq.output ${n_unique_kmers} ${n_tax_reads} resources/pathogensFound.very_inclusive.tab &> logs/FILTER_KRAKENUNIQ_OUTPUT/${sample_name}.log;
        cut -f7 results/KRAKENUNIQ/${sample_name}/krakenuniq.output.pathogens | tail -n +2 > results/KRAKENUNIQ/${sample_name}/taxID.pathogens
done
```

In summary, this rule uses a python script to filter the output from KrakenUniq according to a specified minimum amount of unique kmers and minimum amount of taxReads (reads specific to the taxonomic clade of this species) and extract the information for the species that meet these criteria. 

Please run this code and do not forget to change your account name:

```bash
sbatch KrakenUniq_Filter.sh --account=your_user_account
```

And let's check the output of one sample, `sample1`:

```bash
ls -ltrh results/KRAKENUNIQ/sample1
```

Let's check the output file:

```bash
less results/KRAKENUNIQ/sample1/krakenuniq.output.filtered
```

### Creating an abundance matrix from the KrakenUniq outputs

At last, we will combine the filtered outputs and create an abundance matrix. In this part, the inputs are the filtered krakenuniq outputs of all the samples.

```
rule KrakenUniq_AbundanceMatrix:
    output:
        out_dir=directory("results/KRAKENUNIQ_ABUNDANCE_MATRIX"),
        unique_species="results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt",
        unique_species_names="results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_names_list.txt",
        abundance_matrix="results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix.txt",
        abundance_plot="results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_absolute_abundance_heatmap.pdf",
    input:
        expand("results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered", sample=SAMPLES),
    log:
        "logs/KRAKENUNIQ_ABUNDANCE_MATRIX/KRAKENUNIQ_ABUNDANCE_MATRIX.log",
    params:
        exe=WORKFLOW_DIR / "scripts/krakenuniq_abundance_matrix.R",
        exe_plot=WORKFLOW_DIR / "scripts/plot_krakenuniq_abundance_matrix.R",
        n_unique_kmers=config["n_unique_kmers"],
        n_tax_reads=config["n_tax_reads"],
    message:
        "KrakenUniq_AbundanceMatrix: COMPUTING KRAKENUNIQ MICROBIAL ABUNDANCE MATRIX"
    shell:
        "Rscript {params.exe} results/KRAKENUNIQ {output.out_dir} {params.n_unique_kmers} {params.n_tax_reads} &> {log};"
        "Rscript {params.exe_plot} {output.out_dir} {output.out_dir} &> {log}"
```

This rule combines krakenuniq abundance output files into one file `results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix.txt`.

Then, it creates a heatmap plot using the abundance data: `results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_absolute_abundance_heatmap.pdf`

Please run this script using this piece of code:

```bash
sbatch KrakenUniq_AbundanceMatrix.sh --account=your_account_name

```

Then let's check the output:

```bash
ls results/KRAKENUNIQ_ABUNDANCE_MATRIX/
```
