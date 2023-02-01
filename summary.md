To plot the authentication score and summarise the output of the workflow, we use the following code:


```
rule Plot_Authentication_Score:
"""Plot authentication score"""
    output:
        heatmap="results/overview_heatmap_scores.pdf",    
    input:
        scores=expand("results/AUTHENTICATION/.{sample}_done",sample=SAMPLES)
    message:
        "Plot_Authentication_Score: PLOTTING HEATMAP OF AUTHENTICATION SCORES"
    params:
        exe=WORKFLOW_DIR / "scripts/plot_score.R",
    shell:
        "Rscript {params.exe} results/AUTHENTICATION $(dirname {output.heatmap}) &> {log}"
```

Here is a simplified version of this code:

```bash
"Rscript {params.exe} results/AUTHENTICATION $(dirname {output.heatmap})"

```