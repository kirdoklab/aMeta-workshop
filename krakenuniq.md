# KrakenUniq

## KrakenUniq database

A full NT database for KrakenUniq is [available for download]
(https://www.biorxiv.org/node/2777891.external-links.html) through SciLifeLab. Other databases are available online and for this workshop, we will be using [MinusB](https://benlangmead.github.io/aws-indexes/k2). The database is installed on the cluster.

KrakenUniq is relatively fast and allows to screen aDNA samples against as large database as possible. The KrakenUniq paper suggests that KrakenUniq is not less accurate compared to alignment tools such as BLAST and MEGAN. 

KrakenUniq is very handy as it provides k-mer coverage information that is equivalent to "breadth of coverage" that one can extract via alignment. 

Here we show an example of how to run a sample using the full NT KrakenUniq database:


## Run KrakenUniq



```bash
DBNAME=/proj/mnt/NEOGENE4/share/KrakenUniqDatabase
PATH_TO_KRAKENUNIQ=/usr/local/sw/anaconda3/envs/aMeta/KrakenUniq
time $PATH_TO_KRAKENUNIQ/./krakenuniq --db $DBNAME --fastq-input sample_name.fastq.gz --threads 80 --output sample_name.fastq.gz_sequences.krakenuniq_Full_NT --report-file sample_name.fastq.gz_krakenuniq.output_Full_NT --gzip-compressed --only-classified-out

```

Here we demonstrate how to filter the KrakenUniq output using the kmers=1000 threshold, extract sequence IDs associated with the taxIDs that passed the kmers=1000 threshold and visualize the abundances with Krona:



## Screening workflow


```bash
echo "STEP1: PREPARING PIPELINE"
module load bioinfo-tools; module load perl; module load bowtie2; module load samtools; module load mapDamage/2.0.9; module load gnuparallel; module load Kraken2; module load FastQC
FASTQ_FILE=$1; cd $SNIC_TMP
if [ ! -d $SNIC_TMP/${FASTQ_FILE}_PipelineOutput ];
then
cp /proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/fastq_files/${FASTQ_FILE} $SNIC_TMP; mkdir $SNIC_TMP/${FASTQ_FILE}_PipelineOutput
fi
if [ ! -d $SNIC_TMP/DBDIR ];
then
cp -r /proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/KrakenUniq/krakenuniq/KrakenUniq/DBDIR $SNIC_TMP; cp -r /proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/PathoGenome $SNIC_TMP
fi
SAMPLE=$SNIC_TMP/${FASTQ_FILE}; DBNAME=$SNIC_TMP/DBDIR; OUTPUT=$SNIC_TMP/${FASTQ_FILE}_PipelineOutput; PATH_TO_PATHO_GENOME=$SNIC_TMP/PathoGenome; PATHO_GENOME=$PATH_TO_PATHO_GENOME/library.pathogen.fna
PATH_TO_KRAKENUNIQ=/proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/KrakenUniq/krakenuniq/KrakenUniq; PATH_TO_SCRIPTS=/proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/scripts
KRAKEN2_DB=/proj/snic2018-8-150/uppstore2018095/private/NBIS_Demo/DBDIR_Kraken2_MicrobialNT
 
printf "\n"; echo "STEP2: RUNNING QUALITY CONTROL WITH FASTQC"
if [ ! -f $OUTPUT/*_fastqc.zip ];
then
time fastqc ${SAMPLE} --outdir $OUTPUT
else
echo "SKIPPING RUNNING QUALITY CONTROL WITH FASTQC BECAUSE OUTPUT EXISTS"
fi
 
printf "\n"; echo "STEP3: RUNNING KRAKEN2"
if [ ! -f $OUTPUT/kraken2.output ];
then
time kraken2 ${SAMPLE} --db $KRAKEN2_DB --gzip-compressed --threads 10 --classified-out $OUTPUT/classified_sequences.kraken2 --unclassified-out $OUTPUT/unclassified_sequences.kraken2 --output $OUTPUT/sequences.kraken2 --use-names --report $OUTPUT/kraken2.output
else
echo "SKIPPING RUNNING KRAKEN2 BECAUSE OUTPUT EXISTS"
fi
 
printf "\n"; echo "STEP4: RUNNING KRAKEN-UNIQ"
if [ ! -f $OUTPUT/krakenuniq.output ];
then
time $PATH_TO_KRAKENUNIQ/./krakenuniq --db $DBNAME --fastq-input $SAMPLE --threads 10 --output $OUTPUT/sequences.krakenuniq --report-file $OUTPUT/krakenuniq.output --gzip-compressed --only-classified-out
else
echo "SKIPPING RUNNING KRAKEN-UNIQ BECAUSE OUTPUT EXISTS"
fi
 
printf "\n"; echo "STEP5: FILTERING KRAKEN-UNIQ OUTPUT"
cd $OUTPUT; Rscript $PATH_TO_SCRIPTS/pipeline.R
cut -f7 $OUTPUT/krakenuniq.output.pathogens | tail -n +2 > $OUTPUT/taxID.pathogens
 
printf "\n"; echo "STEP6: COMPUTING MEAN READ LENGTH"
if [ ! -f $OUTPUT/mean.reads.length ];
then
cat $OUTPUT/taxID.pathogens | parallel "${PATH_TO_KRAKENUNIQ}/./krakenuniq-extract-reads {} $OUTPUT/sequences.krakenuniq ${SAMPLE} > $OUTPUT/{}.temp.fq"
echo "MeanReadLength" > $OUTPUT/mean.reads.length; cd $OUTPUT
for i in $(cat taxID.pathogens); do awk '{if(NR%4==2) print length($1)}' ${i}.temp.fq | awk '{ sum += $0 } END { if (NR > 0) print sum / NR }' >> mean.reads.length; done; rm *.temp.fq
else
echo "SKIPPING COMPUTING MEAN READ LENGTH BECAUSE OUTPUT EXISTS"
fi
 
printf "\n"; echo "STEP7: GENERATING FINAL KRAKENUNIQ OUTPUT"
paste krakenuniq.output.pathogens mean.reads.length > krakenuniq.output.pathogens_with_mean_read_length
cat krakenuniq.output.pathogens_with_mean_read_length
 
printf "\n"; echo "STEP8: PERFORMING ALIGNMENT TO PATHO-GENOME"
time bowtie2 --large-index -x $PATHO_GENOME --threads 10 --end-to-end --very-sensitive -U $SAMPLE | samtools view -bS -q 1 -h -@ 10 - > $OUTPUT/test_sample_AlignedToSpecies.bam
samtools sort $OUTPUT/test_sample_AlignedToSpecies.bam -@ 10 > $OUTPUT/test_sample_AlignedToSpecies.sorted.bam; samtools index $OUTPUT/test_sample_AlignedToSpecies.sorted.bam
samtools markdup -r $OUTPUT/test_sample_AlignedToSpecies.sorted.bam $OUTPUT/test_sample_AlignedToSpecies.sorted.dedup.bam; samtools index $OUTPUT/test_sample_AlignedToSpecies.sorted.dedup.bam
 
printf "\n"; echo "STEP9: EXTRACTING ALIGNMENTS FOR CANDIDATES"
cat $OUTPUT/taxID.pathogens | parallel "grep -w {} ${PATH_TO_PATHO_GENOME}/seqid2taxid.pathogen.map | cut -f1 > ${OUTPUT}/{}.seq.ids"
for i in $(cat $OUTPUT/taxID.pathogens); do samtools view -bh $OUTPUT/test_sample_AlignedToSpecies.sorted.dedup.bam -@ 10 $(cat $OUTPUT/${i}.seq.ids | tr "\n" " ") > $OUTPUT/${i}.output.bam; done
 
printf "\n"; echo "STEP10: RUNNING MAPDAMAGE ANCIENT STATUS ANALYSIS"
find . -name '*.output.bam' | parallel "mapDamage -i {} -r ${PATHO_GENOME} --merge-reference-sequences -d ${OUTPUT}/results_{}"
 
printf "\n"; echo "STEP11: ASSIGN ANCIENT STATUS"
Rscript $PATH_TO_SCRIPTS/ancient_status.R 0.05 0.9 $OUTPUT
 
printf "\n"; echo "STEP12: COMPUTE DEPTH AND BREADTH OF COVERAGE FROM ALIGNMENTS"
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


### Kraken2 vs KrakenUniq


Kraken2 is incredibly fast but unfortunately has a high False Positive Rate. KrakenUniq reduces the fraction of false discoveries (compared to Kraken2) by reporting the number of unique k-mers that is analogous to breadth of coverage information. Therefore filtering KrakenUniq output by number of unique k-mers, we get a more reliable list of candidates compared to Kraken2 (the output from Kraken2 can only be filtered by microbial abundance, i.e. depth of coverage). Considering filtered KrakenUniq output as a ground truth, we investigated how many records from Kraken2 output should be retained in order to capture all the species identified with KareknUniq, i.e. how long down the Kraken2 list one needs to go in order to detect all the species reported by KrakenUniq.

What we can conclude after having screened ~100 samples is that a typical sample contains ~100 species reported (after filtering) by KrakenUniq. This corresponds to approximately ~1000 species in Kraken2 output (or 0.01% assigned reads). This implies that if we select top 1000 species in Kraken2, 900 species will probably be false positives and only 10% of all species will be true positives. If we consider KrakenUniq output as a ground truth.
