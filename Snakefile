SAMPLES = ['SRR2090174']

rule all:
    input:
        expand("outputs/megahit/{sample}.contigs.fa", sample = SAMPLES)

rule download_reads:
    output: 
        r1='inputs/raw/{sample}_1.fastq.gz',
        r2='inputs/raw/{sample}_2.fastq.gz'
    shell:'''
    wget -O {output.r1} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/004/SRR2090174/SRR2090174_1.fastq.gz
    wget -O {output.r2} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/004/SRR2090174/SRR2090174_2.fastq.gz    
    '''

rule download_cmag:
    output: "inputs/cmag/GCA_001430905.1_ASM143090v1_genomic.fna.gz"
    shell:'''
    wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/430/905/GCA_001430905.1_ASM143090v1/GCA_001430905.1_ASM143090v1_genomic.fna.gz
    '''

rule adapter_trim:
    input:
        r1 = "inputs/raw/{sample}_1.fastq.gz",
        r2 = 'inputs/raw/{sample}_2.fastq.gz',
        adapters = 'inputs/adapters.fa'
    output:
        r1 = 'outputs/trim/{sample}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{sample}_R2.trim.fq.gz',
        o1 = 'outputs/trim/{sample}_o1.trim.fq.gz',
        o2 = 'outputs/trim/{sample}_o2.trim.fq.gz'
    conda: 'env.yml'
    shell:'''
     trimmomatic PE {input.r1} {input.r2} \
             {output.r1} {output.o1} {output.r2} {output.o2} \
             ILLUMINACLIP:{input.adapters}:2:0:15 MINLEN:25  \
             LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2
    '''

rule kmer_trim_reads:
    input: 
        "outputs/trim/{sample}_R1.trim.fq.gz", 
        "outputs/trim/{sample}_R2.trim.fq.gz"
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'env.yml'
    shell:'''
    interleave-reads.py {input} | 
        trim-low-abund.py -C 3 -Z 18 -M 30e9 -V - -o {output}
    '''

rule sgc_bin_queries:
    input: 
        conf = "inputs/conf/{sample}_conf.yml",
        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz",
        cmag = ["inputs/cmag/GCA_001430905.1_ASM143090v1_genomic.fna.gz"]
    output: "{sample}_k31_r1_search_oh0/GCA_001430905.1_ASM143090v1_genomic.fna.gz.cdbg_ids.reads.fa.gz"
    shell:'''
    python -m spacegraphcats {input.conf} extract_contigs extract_reads --nolock
    '''

rule index_cmag:
    input: "inputs/cmag/GCA_001430905.1_ASM143090v1_genomic.fna.gz"
    output: "inputs/cmag/GCA_001430905.1_ASM143090v1_genomic.fna.gz.bwt"
    conda: 'env.yml'
    shell:'''
    bwa index {input}
    '''
    
rule map_nbhd_reads:
    input: 
        nbhd_reads="{sample}_k31_r1_search_oh0/GCA_001430905.1_ASM143090v1_genomic.fna.gz.cdbg_ids.reads.fa.gz",
        cmag="inputs/cmag/GCA_001430905.1_ASM143090v1_genomic.fna.gz",
        index="inputs/cmag/GCA_001430905.1_ASM143090v1_genomic.fna.gz.bwt"
    output: "outputs/map_nbhd_reads/{sample}_nbhd.sam"
    conda: 'env.yml'
    shell:'''
    bwa mem -t 4 {input.cmag} {input.nbhd_reads} > {output}
    '''

rule sam_unmapped_reads:
    input: "outputs/map_nbhd_reads/{sample}_nbhd.sam"
    output: "outputs/map_nbhd_reads/{sample}_unmapped.sam"
    conda: 'env.yml'
    shell:'''
    samtools view -f 4 {input} > {output}
    ''' 

rule unmapped_reads_to_fastq:
    input: "outputs/map_nbhd_reads/{sample}_unmapped.sam"
    output: "outputs/map_nbhd_reads/{sample}_unmapped.fa"
    conda: 'env.yml'
    shell:'''
    samtools fasta {input} > {output}
    '''

rule megahit_unmapped_reads:
    input: "outputs/map_nbhd_reads/{sample}_unmapped.fa"
    output: "outputs/megahit/{sample}.contigs.fa"
    conda: 'env.yml'
    shell:'''
    megahit -r {input} --min-contig-len 142 \
        --out-dir {wildcards.sample}_megahit \
        --out-prefix {wildcards.sample}
    mv {wildcards.sample}_megahit/{wildcards.sample}.contigs.fa {output}
    rm -rf {wildcards.sample}_megahit
    '''
