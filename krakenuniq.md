
## KrakenUniq database

KrakenUniq is relatively fast and allows to screen aDNA samples against a database as large as possible. The KrakenUniq paper suggests that KrakenUniq is not less accurate compared to alignment tools such as `BLAST` and `MEGAN`.

A full NT database for KrakenUniq is [available for download](https://www.biorxiv.org/node/2777891.external-links.html) through SciLifeLab.

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
    shell:
        "krakenuniq --preload --db {params.DB} --fastq-input {input.fastq} --threads {threads} --output {output.seqs} --report-file {output.report} --gzip-compressed --only-classified-out &> {log}"
```

Here is a simplified version of this code:

```bash
DBNAME=db/
PATH=${PATH}:/truba/home/egitim/miniconda3/envs/aMeta/bin/
krakenuniq --preload --db $DBNAME --fastq-input ${sample_name} --threads 4 --output ${sample_name}.sequences.krakenuniq --report-file ${sample_name}.krakenuniq.output --gzip-compressed --only-classified-out &> logs/KRAKENUNIQ/${sample_name}.log

```

### Filtering the KrakenUniq output

After running KrakenUniq, the next step is to filter the output. The input is the `krakenuniq.output` generated in the previous step. The outputs will be a filtered KrakenUniq abundance file called `krakenuniq.output.filtered`, a pathogens file `krakenuniq.output.pathogens`, and taxids of these pathogens `taxID.pathogens`.

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
