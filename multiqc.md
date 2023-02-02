MultiQC is a useful tool for getting a more extensive quality report.

```
rule MultiQC:
    """Run MultiQC"""
    output:
        html="results/MULTIQC/multiqc_report.html",
    input:
        unpack(multiqc_input),
    log:
        "logs/MULTIQC/MULTIQC.log",
    conda:
        "../envs/multiqc.yaml"
    params:
        config=os.path.join(WORKFLOW_DIR, "envs", "multiqc_config.yaml"),
    envmodules:
        *config["envmodules"]["multiqc"],
    benchmark:
        "benchmarks/MULTIQC/MULTIQC.benchmark.txt"
    message:
        "MultiQC: COMBINING QUALITY CONTROL METRICS WITH MULTIQC"
    shell:
        'echo {input} | tr " " "\n" > {output.html}.fof;'
        "multiqc -c {params.config} -l {output.html}.fof --verbose --force --outdir results/MULTIQC &> {log}"
```

