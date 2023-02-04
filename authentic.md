---
title: "Authentication of the pathogens"

---

In this part, we will start authenticating the pathogens found by the `krakenuniq` tool. For this example part, we will focus on the pathogens particular for `sample1`.

The main logic in this section, is to extract DNA reads aligned to one specific pathogen, and run authentication commands.

In the `krakenuniq` part, we created a file called `taxID.pathogens`.

Let's check this file:

```bash
results/KRAKENUNIQ/sample1/taxID.pathogens
```

```
mkdir -p results/AUTHENTICATION/sample1/13373/
```

Afterwards, we will extract the node name from the `krakenuniq` database. Let's check the output:

```bash

less results/AUTHENTICATION/sample1/13373/node_list.txt

```

THis pathogen name is Burkholderia mallei
