SAMPLES, = glob_wildcards("data/samples/{sample}.1.fastq.gz")

include: "rules/multiqc.smk"
include: "rules/map.smk"

rule all:
    input: [expand("data/salmon/{sample}/quant.sf", sample=SAMPLES), "data/multiqc/multiqc_report.html"]
