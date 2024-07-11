SAMPLES, = glob_wildcards("data/samples/{sample}.1.fastq.gz")

include: "rules/multiqc.smk"
include: "rules/map_salmon.smk"
include: "rules/map_star.smk"

rule all:
    input: [expand("data/salmon/{sample}/quant.sf", sample=SAMPLES),expand("data/star/{sample}_Aligned.sortedByCoord.out.bam",sample=SAMPLES), "data/multiqc/multiqc_report.html"]
