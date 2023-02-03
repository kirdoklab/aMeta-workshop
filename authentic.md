---
title: "Authentication steps"
---

## Create a directory per TaxID per sample

This rule is actually a snakemake checkpoint and you can learn more about it in the section "How to understand snakemake?". Basically, it creates a directory per TaxID found after filtering per sample.

```
checkpoint Create_Sample_TaxID_Directories:
    """Create taxid directory
    For a sample, create taxid for each entry in krakenuniq output
    taxID.pathogens. Downstream rules use the taxid directories as
    input, but it is not known beforehand which these are; they are
    determined by the finds in krakenuniq.
    """
    input:
        pathogens="results/KRAKENUNIQ/{sample}/taxID.pathogens",
    output:
        done="results/AUTHENTICATION/{sample}/.extract_taxids_done",
    log:
        "logs/CREATE_SAMPLE_TAXID_DIRECTORIES/{sample}.log",
    params:
        dir=lambda wildcards: f"results/AUTHENTICATION/{wildcards.sample}",
    shell:
        "mkdir -p {params.dir}; "
        "while read taxid; do mkdir -p {params.dir}/$taxid; touch {params.dir}/$taxid/.done; done<{input.pathogens};"
        "touch {output.done}"
```

```bash
for fastq_file in results/CUTADAPT_ADAPTER_TRIMMING/*; do
  sample=$(basename "${fastq_file}" .fastq.gz)
  mkdir -p results/AUTHENTICATION/${sample} logs/AUTHENTICATION;
  while read taxid; do mkdir -p results/AUTHENTICATION/${sample}/${taxid}; touch results/AUTHENTICATION/${sample}/${taxid}/.done; done<results/KRAKENUNIQ/${sample}/taxID.pathogens;
  touch results/AUTHENTICATION/${sample}/.extract_taxids_done
done
```

In summary, this command loops through all the samples and creates an authentication folder for each of them. Then it reads the taxID.pathogens file containing the list of pathogens that were found and creates a folder for each taxID per sample.

## Make_Node_list 

This rule creates a file called node_list.txt containing the name of the species for that taxID.

WARNING: We won't be able to run this rule for the course because we don't have access to the KrakenUniq database for this workshop.

```
rule Make_Node_List:
    """Generate a list of species names for a taxonomic identifier"""
    input:
        dirdone="results/AUTHENTICATION/{sample}/{taxid}/.done",
    output:
        node_list="results/AUTHENTICATION/{sample}/{taxid}/node_list.txt",
    params:
        tax_db=config["krakenuniq_db"],
    log:
        "logs/MAKE_NODE_LIST/{sample}_{taxid}.log",
    shell:
        "awk -v var={wildcards.taxid} '{{ if($1==var) print $0 }}' {params.tax_db}/taxDB | cut -f3 > {output.node_list}"
```

Here is a shell version of this rule:

```bash
tax_db=resources/DBDIR_KrakenUniq_MicrobialNT
for fastq_file in results/CUTADAPT_ADAPTER_TRIMMING/*; do
  sample=$(basename "${fastq_file}" .fastq.gz)
  while read taxid; do
    awk -v var=$taxid '{{ if($1==var) print $0 }}' $tax_db/taxDB | cut -f3 > results/AUTHENTICATION/$sample/$taxid/node_list.txt
  done < results/KRAKENUNIQ/$sample/taxID.pathogens
done
```

In summary, this rule looks for each sample and each taxid into the KrakenUniq database to find the corresponding latin name for each taxID and write this lating name into the file node_list.txt.
