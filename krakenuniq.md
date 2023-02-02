
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


### Filtering the output

After running KrakenUniq, the next step is to filter the output. The input is the `krakenuniq.output` generated in the previous step.

```
rule Filter_KrakenUniq_Output:
    output:
        filtered="results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered",
        pathogens="results/KRAKENUNIQ/{sample}/krakenuniq.output.pathogens",
        pathogen_tax_id="results/KRAKENUNIQ/{sample}/taxID.pathogens",
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

Breadth and depth of coverage filters 
default thresholds are very conservative, can be tuned by users

```bash
n_unique_kmers: 1000
n_tax_reads: 200

cd $OUTPUT; Rscript $PATH_TO_SCRIPTS/pipeline.R
cut -f7 $OUTPUT/krakenuniq.output.pathogens | tail -n +2 > $OUTPUT/taxID.pathogens
```

```
cat $OUTPUT/taxID.pathogens | parallel "${PATH_TO_KRAKENUNIQ}/./krakenuniq-extract-reads {} $OUTPUT/sequences.krakenuniq ${SAMPLE} > $OUTPUT/{}.temp.fq"
echo "MeanReadLength" > $OUTPUT/mean.reads.length; cd $OUTPUT
for i in $(cat taxID.pathogens); do awk '{if(NR%4==2) print length($1)}' ${i}.temp.fq | awk '{ sum += $0 } END { if (NR > 0) print sum / NR }' >> mean.reads.length; done; rm *.temp.fq

```

```bash
paste krakenuniq.output.pathogens mean.reads.length > krakenuniq.output.pathogens_with_mean_read_length
cat krakenuniq.output.pathogens_with_mean_read_length
```


### KrakenUniq to Krona

We then visualize our filtered output using Krona:

```
rule KrakenUniq2Krona:
    output:
        tax_ids="results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered_taxIDs_kmers1000.txt",
        seqs="results/KRAKENUNIQ/{sample}/sequences.krakenuniq_kmers1000.txt",
        krona="results/KRAKENUNIQ/{sample}/sequences.krakenuniq_kmers1000.krona",
        html="results/KRAKENUNIQ/{sample}/taxonomy.krona.html",
    input:
        report="results/KRAKENUNIQ/{sample}/krakenuniq.output.filtered",
        seqs="results/KRAKENUNIQ/{sample}/sequences.krakenuniq",
    log:
        "logs/KRAKENUNIQ2KRONA/{sample}.log",
    conda:
        "../envs/krona.yaml"
    envmodules:
        *config["envmodules"]["krona"],
    params:
        exe=WORKFLOW_DIR / "scripts/krakenuniq2krona.py",
        DB=f"--tax {config['krona_db']}" if "krona_db" in config else "",
    benchmark:
        "benchmarks/KRAKENUNIQ2KRONA/{sample}.benchmark.txt"
    message:
        "KrakenUniq2Krona: VISUALIZING KRAKENUNIQ RESULTS WITH KRONA FOR SAMPLE {input.report}"
    shell:
        "{params.exe} {input.report} {input.seqs} &> {log}; "
        "cat {output.seqs} | cut -f 2,3 > {output.krona}; "
        "ktImportTaxonomy {output.krona} -o {output.html} {params.DB} &>> {log}"
```

### AbundanceMatrix

Lastly, we create an Abundance matrix:

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



```bash
printf "\n"; echo "PERFORMING ALIGNMENT TO PATHO-GENOME"
time bowtie2 --large-index -x $PATHO_GENOME --threads 10 --end-to-end --very-sensitive -U $SAMPLE | samtools view -bS -q 1 -h -@ 10 - > $OUTPUT/test_sample_AlignedToSpecies.bam
samtools sort $OUTPUT/test_sample_AlignedToSpecies.bam -@ 10 > $OUTPUT/test_sample_AlignedToSpecies.sorted.bam; samtools index $OUTPUT/test_sample_AlignedToSpecies.sorted.bam
samtools markdup -r $OUTPUT/test_sample_AlignedToSpecies.sorted.bam $OUTPUT/test_sample_AlignedToSpecies.sorted.dedup.bam; samtools index $OUTPUT/test_sample_AlignedToSpecies.sorted.dedup.bam
 
printf "\n"; echo "EXTRACTING ALIGNMENTS FOR CANDIDATES"
cat $OUTPUT/taxID.pathogens | parallel "grep -w {} ${PATH_TO_PATHO_GENOME}/seqid2taxid.pathogen.map | cut -f1 > ${OUTPUT}/{}.seq.ids"
for i in $(cat $OUTPUT/taxID.pathogens); do samtools view -bh $OUTPUT/test_sample_AlignedToSpecies.sorted.dedup.bam -@ 10 $(cat $OUTPUT/${i}.seq.ids | tr "\n" " ") > $OUTPUT/${i}.output.bam; done
 
printf "\n"; echo "RUNNING MAPDAMAGE ANCIENT STATUS ANALYSIS"
find . -name '*.output.bam' | parallel "mapDamage -i {} -r ${PATHO_GENOME} --merge-reference-sequences -d ${OUTPUT}/results_{}"
 
printf "\n"; echo "ASSIGN ANCIENT STATUS"
Rscript $PATH_TO_SCRIPTS/ancient_status.R 0.05 0.9 $OUTPUT

```


```bash
printf "\n"; echo "COMPUTE DEPTH AND BREADTH OF COVERAGE FROM ALIGNMENTS"
echo "NUMBER_OF_READS" > DepthOfCoverage.${FASTQ_FILE}.txt; echo "GENOME_LENGTH" > GenomeLength.${FASTQ_FILE}.txt; echo "BREADTH_OF_COVERAGE" > BreadthOfCoverage.${FASTQ_FILE}.txt
for j in $(cut -f7 final_output.txt | sed '1d')
do
echo "Organism $j"
if [ -s ${j}.output.bam ] && [ "$(samtools view ${j}.output.bam | wc -l)" -ne "0" ];
then
samtools sort ${j}.output.bam > ${j}.output.sorted.bam; samtools depth ${j}.output.sorted.bam | cut -f1 | uniq > Genomes_${j}.txt
GENOME_LENGTH=$(grep -wFf Genomes_${j}.txt $PATH_TO_PATHO_GENOME/GenomeLength.txt | cut -f2 | awk '{ sum += $1; } END { print sum; }')
if [ -s ${j}.output.sorted.bam ];
then
NUMBER_OF_READS=$(samtools view ${j}.output.sorted.bam | wc -l); NUMBER_OF_COVERED_POSITIONS=$(samtools depth ${j}.output.sorted.bam | wc -l)
else
echo "NA" >> DepthOfCoverage.${FASTQ_FILE}.txt; echo "NA" >> BreadthOfCoverage.${FASTQ_FILE}.txt; echo "NA" >> GenomeLength.${FASTQ_FILE}.txt
continue
fi
echo $NUMBER_OF_READS >> DepthOfCoverage.${FASTQ_FILE}.txt; echo $GENOME_LENGTH >> GenomeLength.${FASTQ_FILE}.txt
echo "scale=10 ; ($NUMBER_OF_COVERED_POSITIONS / $GENOME_LENGTH)*100" | bc >> BreadthOfCoverage.${FASTQ_FILE}.txt
else
echo "NA" >> DepthOfCoverage.${FASTQ_FILE}.txt; echo "NA" >> BreadthOfCoverage.${FASTQ_FILE}.txt; echo "NA" >> GenomeLength.${FASTQ_FILE}.txt
continue
fi
done; rm Genomes_*.txt; paste final_output.txt DepthOfCoverage.${FASTQ_FILE}.txt GenomeLength.${FASTQ_FILE}.txt BreadthOfCoverage.${FASTQ_FILE}.txt > final_output_corrected.txt
 
cp -r $OUTPUT /proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo; cat final_output_corrected.txt
printf "\n"; echo "PIPELINE FINISHED SUCCESSFULLY"

```

