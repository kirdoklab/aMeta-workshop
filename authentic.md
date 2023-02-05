---
title: "Authentication of the pathogens"

---

In this part, we will start authenticating the pathogens found by the `krakenuniq` tool. This part is complex and needs database dependencies to work. So in this example part, we will focus on the pathogen (taxid 13373) particular for `sample1`.

Before working on this part, lets export the `PATH` variable to check the output files:

```bash
export PATH=${PATH}:/truba/home/egitim/miniconda3/envs/aMeta/bin/
```

The main logic in this section, is to extract DNA reads aligned to one specific pathogen, and run authentication commands.

Let's first run our script, and check what is going in here.

```bash
sbatch Authentic.sh --account=egitim  
```

Let's go step by step.

In the `krakenuniq` part, we created a file called `taxID.pathogens`.

Let's check this file:

```bash
less /truba/home/egitim/aMeta/results/KRAKENUNIQ/sample1/taxID.pathogens
```

Afterwards, we will extract the node name from the `krakenuniq` database. We can not show the output, because it needs the big krakenuniq database.

Let's check the output:

```bash

less /truba/home/egitim/aMeta/results/AUTHENTICATION/sample1/13373/node_list.txt

```

THis pathogen name is *Burkholderia mallei*. Over the next steps, we will extract DNA reads assigned to this pathogen, and we will create authenticity metrics.

Then we will use `MaltExtract` and `postprocessing.AMPS.r` tools to extract DNA reads assigned to this pathogen, from the rma6 file of the sample1.

Let's check this folder:

```bash

ls /truba/home/egitim/aMeta/results/AUTHENTICtput//sample1/13373/MaltExtract_out/ancient
```

Then we will extract the sequence name of the reference sequence of the bacteria from the database:

```bash
/truba/home/egitim/aMeta/results/AUTHENTICATION/sample1/13373/name_list.txt
```

Then we extract alignment entries from the malt `sam` file using this sequence ID. Let's check the output file:

```bash

```