SAMPLES, = glob_wildcards("data/samples/{sample}.1.fastq.gz")

include: "rules/multiqc.smk"
include: "rules/map_salmon.smk"
include: "rules/map_star.smk"

rule all:
    input: [expand("data/salmon/{sample}/quant.sf", sample=SAMPLES),expand("data/star/{sample}_deduplicated.bam",sample=SAMPLES), "data/multiqc/multiqc_report.html",expand("data/tracks/{sample}_fwd.bw",sample=SAMPLES),expand("data/tracks/{sample}_rev.bw",sample=SAMPLES)]
