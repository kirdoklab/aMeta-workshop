

For a sample, create taxid for each entry in krakenuniq output taxID.pathogens. Downstream rules use the taxid directories as input, but it is not known beforehand which these are; they are determined by the finds in krakenuniq.


```
rule Authentication:
    """Run Authentication"""
    output:
        done="results/AUTHENTICATION/{sample}/.extract_taxids_done",
    input:
        pathogens="results/KRAKENUNIQ/{sample}/taxID.pathogens",
    shell:
        "mkdir -p {params.dir}; "
        "while read taxid; do mkdir -p {params.dir}/$taxid; touch {params.dir}/$taxid/.done; done<{input.pathogens};"
        "touch {output.done}"
```

Here is a simplified version of this code:

```bash
"mkdir -p {params.dir}; "
        "while read taxid; do mkdir -p {params.dir}/$taxid; touch {params.dir}/$taxid/.done; done<{input.pathogens};"
        "touch {output.done}"
```
        
        
 ```
for BASE in ${SAMPLES}
do
	sbatch --time=12:00:00 --job-name=auth_${BASE} --ntasks-per-node=${THREADS} -A ${PROJECT} --mail-user=${MAIL} --mail-type=${NOTIFICATION} ${MALTOUTPUT}/to_run_authenticate.sh ${BASE} ${MALTOUTPUT}/${BASE} ${FILE}
done

