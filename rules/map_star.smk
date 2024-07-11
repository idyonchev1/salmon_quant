rule map_star:
        input:["data/trimmed/trimmed_{sample}.1.fastq","data/trimmed/trimmed_{sample}.2.fastq"]
        output:"data/star/{sample}_Aligned.sortedByCoord.out.bam"
        log:"logs/star/{sample}.log"
        threads:32
        resources:
            mem_mb=50000,
            time="48:00:00"
        shell:
            """
            ulimit -n 65535
            STAR --runThreadN {threads} --genomeDir star_index --readFilesIn {input[0]} {input[1]} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "data/star/{wildcards.sample}_"
            samtools index {output}
            """
