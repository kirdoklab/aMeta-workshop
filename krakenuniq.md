
## KrakenUniq database

KrakenUniq is relatively fast and allows to screen aDNA samples against as large database as possible. The KrakenUniq paper suggests that KrakenUniq is not less accurate compared to alignment tools such as `BLAST` and `MEGAN`.

A full NT database for KrakenUniq is [available for download](https://www.biorxiv.org/node/2777891.external-links.html) through SciLifeLab. Other databases are available online and for this workshop, we use [Standard-8](https://benlangmead.github.io/aws-indexes/k2). The database is installed on the cluster.

## KrakenUniq workflow

### Run KrakenUniq

**Please note that it will not be possible to run this first step during this workshop due to issues with storage space. The output files will be provided by the course leaders.**

To start, we need to set up the path to the database– we call it `DBNAME`.

This is the rule in the workflow/rules/krakenuniq.smk file. Input files are retrieved from the output of the Quality Control step.

```
rule KrakenUniq:
    """Run KrakenUniq on trimmed fastq data"""
    output:
        report="results/KRAKENUNIQ/{sample}/krakenuniq.output",
        seqs="results/KRAKENUNIQ/{sample}/sequences.krakenuniq",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
    params:
        DB=config["krakenuniq_db"],
    threads: 10
    log:
        "logs/KRAKENUNIQ/{sample}.log",
    conda:
        "../envs/krakenuniq.yaml"
    envmodules:
        *config["envmodules"]["krakenuniq"],
    benchmark:
        "benchmarks/KRAKENUNIQ/{sample}.benchmark.txt"
    message:
        "KrakenUniq: PERFORMING TAXONOMIC CLASSIFICATION OF SAMPLE {input.fastq} WITH KRAKENUNIQ"
    shell:
        "krakenuniq --preload --db {params.DB} --fastq-input {input.fastq} --threads {threads} --output {output.seqs} --report-file {output.report} --gzip-compressed --only-classified-out &> {log}"
```

Here is a simplified version of this code:

```bash
DBNAME=db/
PATH=${PATH}:/truba/home/egitim/miniconda3/envs/aMeta/bin/
krakenuniq --preload --db $DBNAME --fastq-input ${sample_name} --threads 4 --output ${sample_name}.sequences.krakenuniq --report-file ${sample_name}.krakenuniq.output --gzip-compressed --only-classified-out &> logs/KRAKENUNIQ/${sample_name}.log

```

### Filtering the krakenuniq output

After running KrakenUniq, the next step is to filter the output. The input is the `krakenuniq.output` generated in the previous step. The outputs will be a filtered kraken uniq abundance file `krakenuniq.output.filtered`, pathogens file `krakenuniq.output.pathogens`, and taxids of these pathogens `taxID.pathogens`.

To filter the krakenuniq output, we will be using the breadth and depth of coverage filters. These filters are defined in the `config.yaml` file. Breadth and depth of coverage filters  default thresholds are very conservative, can be tuned by users. 

The rule is this:

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
    conda:
        "../envs/krakenuniq.yaml"
    envmodules:
        *config["envmodules"]["krakenuniq"],
    benchmark:
        "benchmarks/FILTER_KRAKENUNIQ_OUTPUT/{sample}.benchmark.txt"
    message:
        "Filter_KrakenUniq_Output: APPLYING DEPTH AND BREADTH OF COVERAGE FILTERS TO KRAKENUNIQ OUTPUT FOR SAMPLE {input}"
    shell:
        """{params.exe} {input.krakenuniq} {params.n_unique_kmers} {params.n_tax_reads} {input.pathogenomesFound} &> {log}; """
        """cut -f7 {output.pathogens} | tail -n +2 > {output.pathogen_tax_id}"""
        
```

Here is a simplified version of this code:

```bash
n_unique_kmers: 1000
n_tax_reads: 200

cd $OUTPUT; Rscript $PATH_TO_SCRIPTS/pipeline.R
cut -f7 $OUTPUT/krakenuniq.output.pathogens | tail -n +2 > $OUTPUT/taxID.pathogens
cat $OUTPUT/taxID.pathogens | parallel "${PATH_TO_KRAKENUNIQ}/./krakenuniq-extract-reads {} $OUTPUT/sequences.krakenuniq ${SAMPLE} > $OUTPUT/{}.temp.fq"
echo "MeanReadLength" > $OUTPUT/mean.reads.length; cd $OUTPUT
for i in $(cat taxID.pathogens); do awk '{if(NR%4==2) print length($1)}' ${i}.temp.fq | awk '{ sum += $0 } END { if (NR > 0) print sum / NR }' >> mean.reads.length; done; rm *.temp.fq

```

```bash
paste krakenuniq.output.pathogens mean.reads.length > krakenuniq.output.pathogens_with_mean_read_length
cat krakenuniq.output.pathogens_with_mean_read_length
```

Please run this code and do not forget to change your account name:

```bash
sbatch KrakenUniq_Filter.sh --account=egitim
```

And let's check the output of one sample, `sample1`:

```bash
ls -ltrh results/KRAKENUNIQ/sample1
```

Let's check the output file:

```bash
less results/KRAKENUNIQ/sample1/krakenuniq.output.filtered
```

### Creating Abundance Matrix from Krakenuniq outputs

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
    conda:
        "../envs/r.yaml"
    envmodules:
        *config["envmodules"]["r"],
    benchmark:
        "benchmarks/KRAKENUNIQ_ABUNDANCE_MATRIX/KRAKENUNIQ_ABUNDANCE_MATRIX.benchmark.txt"
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
sbatch KrakenUniq_AbundanceMatrix.sh --account=egitim

```

Then let's check the output:

```bash
ls results/KRAKENUNIQ_ABUNDANCE_MATRIX/
```