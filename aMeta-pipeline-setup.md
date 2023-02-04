---
title: "aMeta pipeline setup"
---

## Introduction

Okay, that's all well and good, but when are we really going to learn how to use the snakemake aMeta pipeline?

Now that you know exactly what the different rules of aMeta do, this pipeline is no longer a black box for you and you can start using the official snakemake version. 

But to do so, you should not rush into it but think about different criteria for its installation. 

## Decisions

There are several decisions to be made:

+ Which databases do I need and where will I store them?
+ Where do I install the packages for the conda environment?
+ How do I configure snakemake to use the system on my high-performance computer?

These are some of the questions we will try to answer here. 

### The necessary databases

Depending on the type of project you want to do, there are different databases to download that we have made available to you. 

If your project is only focused on prokaryotes, the KrakenUniq and Bowtie2 microbial databases will be sufficient, but you will still need to download the Bowtie2_FULL_NT to get the files needed to help build the MALT database.

If your project is interested in all organisms, you will need to download the KrakenUniq_FULL_NT and the Bowtie2_FULL_NT. 

Links to the databases and a process for downloading them are available on the README of the [official aMeta GitHub](https://github.com/NBISweden/aMeta). 

If several people in your cluster intend to use the same databases, it is better to store them in a folder accessible to all. It is also important to note that these databases are very large and will require a lot of storage space. 

### Conda packages location

It may be worth thinking about the location of packages and environments installed by conda in several situations: for example, when several users will be using these packages or when the number of files or the memory of your home folder is limited. 

In these cases, you may want to specify the desired location of the conda environment when creating it with the option --prefix.

```bash
conda env create --prefix /path/to/your/conda-dirs workflow/envs/environment.yaml
```

If you do use the option --prefix to choose a location for your conda packages and environments, don't forget to create a hidden ".condarc" file in your home directory containing information about the path like this:

```bash 
pkgs_dirs:
    - /path/to/your/conda-dirs/pkgs
envs_dirs:
    - /path/to/your/conda-dirs/envs
```

### Make snakemake use the queue system of your high-performance computer (HPC, example here with slurm)

Once you have followed the installation steps described in the README file of the [official aMeta GitHub](https://github.com/NBISweden/aMeta), you will need to make your snakemake setup compatible with the queuing system of your high-performance computer (HPC system). Do this step before running the first snakemake command for a project involving real data. If you don't, and run the snakemake pipeline for real data (not the dummy data from the .test folder) you risk to receive complaints from the server maintenance team for overusing the logging node. 

There is a small sentence in the README of the [official aMeta GitHub](https://github.com/NBISweden/aMeta) about advanced profiles for HPC systems:
"For more advanced profiles for different hpc systems, see [Snakemake-Profiles github page](https://github.com/snakemake-profiles)". Click on this link. This will bring you to another GitHub with snakemake profiles available for different cluster environments. 

Here, I will explain how to setup a slurm snakemake profile because it is the most common cluster queuing system I have encountered so far. If your server uses sbatch, squeue and such, it means that it is using a slurm system too. 

In that GitHub page click on the ["slurm" repository](https://github.com/Snakemake-Profiles/slurm). To install that cluster repository, you will need to use cookiecutter from Python. Load the python module of your server and follow the instruction from the Snakemake-Profiles/slurm webpage for a quickstart. 

WARNING: This will prompt a lot of questions about your computing system, so be ready with the user account name, the partition name used as default in your server, etc. Also, think through where you want to install this snakemake profile. Your home directory as suggested in the quickstart ? A private folder ? A shared folder ? It will be useful to fine-tune the config.yaml file for each rule in a different way over time (runtime, memory, partition) and to share it with colleagues so that they don't have to redo the fine-tuning again. But keep in mind that if they modify the config.yaml file you are using while you are running an analysis, it may crash. 

So, we leave you with these thoughts and decisions to make and we think you are good to go! If you need more information, or if something is not clear, don't hesitate to contact us. We aim to improve this website over time. Good luck!
