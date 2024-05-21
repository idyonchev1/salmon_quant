rule cutadapt_trim:
    input:
        ["data/samples/{sample}.1.fastq.gz", "data/samples/{sample}.2.fastq.gz"],
    output:
        fastq1="data/trimmed/trimmed_{sample}.1.fastq",
        fastq2="data/trimmed/trimmed_{sample}.2.fastq",
        qc="data/trimmed/trimmed_{sample}.qc.txt",
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-a AGATCGGAAGAGCA -A AGATCGGAAGAGCA",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--minimum-length 20 -q 30",
    log:
        "logs/cutadapt/{sample}.log",
    threads: 1  # set desired number of threads here
    resources:
        mem_mb=8000,
        time="24:00:00"
    wrapper:
        "v1.32.1/bio/cutadapt/pe"

rule fastqc:
    input:
        ["data/trimmed/trimmed_{sample}.1.fastq", "data/trimmed/trimmed_{sample}.2.fastq"]
    output:
        ["data/fastqc/trimmed_{sample}.1_fastqc.html","data/fastqc/trimmed_{sample}.2_fastqc.html"]
#    conda:
#        "../envs/qc.yml"
    log:
        "logs/fastqc/{sample}.log"
    params:
        threads = "1"
    resources:
        mem_mb=4000,
        time="24:00:00"
    shell:
        "fastqc -o data/fastqc/ -t 1 {input} &> {log}"

rule fastq_screen:
    input:
        ["data/trimmed/trimmed_{sample}.1.fastq","data/trimmed/trimmed_{sample}.2.fastq"]
    output:
        txt="data/fastqc/trimmed_{sample}.fastq_screen.txt",
        png="data/fastqc/trimmed_{sample}.fastq_screen.png"
    params:
        fastq_screen_config="config/fastq_screen.conf",
        subset=100000,
        aligner='bowtie2'
    threads: 8
    resources:
        mem_mb=8000,
        time="24:00:00"
    wrapper:
        "v2.0.0/bio/fastq_screen"

rule multiqc:
    input:
        expand("data/fastqc/trimmed_{sample}.1_fastqc.html", sample=SAMPLES),
        expand("data/fastqc/trimmed_{sample}.2_fastqc.html", sample=SAMPLES),
        expand("data/fastqc/trimmed_{sample}.fastq_screen.txt", sample=SAMPLES),
    output:
        "data/multiqc/multiqc_report.html"
#    conda:
#        "../envs/qc.yml"
    resources:
        mem_mb=2000,
        time="24:00:00"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        "multiqc data/fastqc -o data/multiqc &> {log}"
