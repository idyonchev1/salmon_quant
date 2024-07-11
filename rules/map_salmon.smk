rule salmon_quant_reads:
    input:["data/trimmed/trimmed_{sample}.1.fastq", "data/trimmed/trimmed_{sample}.2.fastq"]
    output:"data/salmon/{sample}/quant.sf"
    log:
        "logs/salmon/{sample}.log"
    threads: 12
    conda:
        "salmon"
    resources:
        mem_mb=4000,
        time="24:00:00"
    shell:
        "salmon quant -i salmon_index/transcripts_index --discardOrphansQuasi --hardFilter --mimicStrictBT2 -l A -1 {input[0]} -2 {input[1]} -p 12 -o data/salmon/{wildcards.sample} &> {log}"

