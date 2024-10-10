rule map_star:
        input:["data/trimmed/trimmed_{sample}.1.fastq","data/trimmed/trimmed_{sample}.2.fastq"]
        output:
            bam=temp("data/star/{sample}_Aligned.sortedByCoord.out.bam"),
            bai=temp("data/star/{sample}_Aligned.sortedByCoord.out.bam.bai")
        log:"logs/star/{sample}.log"
        threads:32
        resources:
            mem_mb=50000,
            time="48:00:00"
        shell:
            """
            ulimit -n 65535
            STAR --runThreadN {threads} --genomeDir star_index --readFilesIn {input[0]} {input[1]} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "data/star/{wildcards.sample}_" &> {log}
            samtools index {output.bam} 2>> {log}
            """

rule deduplicate:
        input:"data/star/{sample}_Aligned.sortedByCoord.out.bam"
        output:"data/star/{sample}_deduplicated.bam"
        log:"logs/star/{sample}_dedup.log"
        threads:4
        resources:
            mem_mb=8000,
            time="24:00:00"
        shell:
            """
            samtools view -b {input} -f 3 -F 1024 -q 30 > {output} 2> {log}
            samtools index {output} 2>>{log}
            """

rule export_stranded_bigwigs:
    input:
        "data/star/{sample}_deduplicated.bam"
    output:
        fwdbam1=temp("data/tracks/{sample}_fwd1.bam"),
        fwdbam2=temp("data/tracks/{sample}_fwd2.bam"),
        fwdbam=temp("data/tracks/{sample}_fwd.bam"),
        revbam1=temp("data/tracks/{sample}_rev1.bam"),
        revbam2=temp("data/tracks/{sample}_rev2.bam"),
        revbam=temp("data/tracks/{sample}_rev.bam"),
        fwdbw="data/tracks/{sample}_fwd.bw",
        revbw="data/tracks/{sample}_rev.bw"
    log:
        "logs/deeptools_{sample}.log"
    threads: 4
    resources:
        mem_mb=10000,
        time="24:00:00"
    conda:
        "cgat-apps"
    shell:
        """
        samtools view -b -f 128 -F 16 {input} > {output.fwdbam1}
        samtools view -b -f 80 {input} > {output.fwdbam2}
        samtools merge -f {output.fwdbam} {output.fwdbam1} {output.fwdbam2}
        samtools index {output.fwdbam}
        samtools view -b -f 144 {input} > {output.revbam1}
        samtools view -b -f 64 -F 16 {input} > {output.revbam2}
        samtools merge -f {output.revbam} {output.revbam1} {output.revbam2}
        samtools index {output.revbam}
        bamCoverage --bam {output.fwdbam} -o {output.fwdbw} -p 4 --normalizeUsing RPKM --exactScaling --binSize 1 --effectiveGenomeSize 2913022398 &> {log} 
        bamCoverage --bam {output.revbam} -o {output.revbw} -p 4 --normalizeUsing RPKM --exactScaling --binSize 1 --effectiveGenomeSize 2913022398 &>> {log} 
        rm data/tracks/{wildcards.sample}*.bai
        """
