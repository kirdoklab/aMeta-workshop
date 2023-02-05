---
title: "Authentication of the pathogens"

---

In this part, we will start authenticating the pathogens found by the `krakenuniq` tool. This part is complex and needs database dependencies to work. So in this example part, we will focus on the pathogen (taxid 13373) particular for `sample1`.

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

Then we will use `MaltExtract` tool to extract DNA reads assigned to this pathogen, from the rma6 file of the sample1.

Let's check this folder:

```bash

ls /truba/home/egitim/aMeta/results/AUTHENTICtput//sample1/13373/MaltExtract_out
```