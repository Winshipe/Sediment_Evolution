import sys
import json
import os

#note that here "depth" ALWAYS refers to sample depth (Cm Below Sea Floor) rather than
#depth of coverage!  I (and instrain) use "coverage" when talking about depth of coverage 

#yes I could put these as config variables, yes that would be better practice
#but no because that takes longer and makes calling the program less "sure"
mag_origin = "M1_3" #"M5_25" #"M5_3"
my_site = "M1" #"M5"
mag_numbers = [216,182,23,34,148,61,88,93,156,191] #[13, 108, 117, 21, 68, 31, 17, 82] #[105,122,180,18,202,248,252,39,3,64,68]
my_depths = [1,3,5,7,10,25]

rule all:
    input:
        expand(
            ["refs/{mag_origin}.{mag}.filtered.fa",
             "trimmed/trimmed_{site}_{depth}_R1.fastq",
            "trimmed/trimmed_{site}_{depth}_R2.fastq",
            "combined_reference_{mag_origin}.fa",
            "combined_reference_{mag_origin}.stb",
            "combined_reference_{mag_origin}.fna",
            "annotations/{mag_origin}.{mag}_eggnog.emapper.annotations",
            "regions_{mag_origin}.tsv",
            "bams/normalized_filtered_{mag_origin}_{depth}.bam",
            "bams/processed_{mag_origin}_{depth}.bam",
            "IS_proc_{mag_origin}_{depth}/output/IS_proc_{mag_origin}_{depth}_genome_info.tsv"],
            site = [my_site],
            mag_origin = [mag_origin],
            depth = my_depths,       
            mag = mag_numbers,
            allow_missing = True
        )


rule combine_references:
    input:
        refs = expand(["refs/{s}.{m}.filtered.fa"], s=mag_origin, m = mag_numbers)
    output:
        comb_fa = "combined_reference_{mag_origin}.fa",
        stb = "combined_reference_{mag_origin}.stb"
    shell: """
        for f in {input.refs}
        do 
        grep '>' $f | awk ' {{ gsub(">","") ; print $1 "\t" "'$f'" }} ' >> {output.stb}
        cat $f >> {output.comb_fa}
        done
    """

rule trim:
    input:
        R1 = "raw_reads/{site}_{depth}_R1.fastq.gz",
        R2 = "raw_reads/{site}_{depth}_R2.fastq.gz",
    output:
        R1 = "trimmed/trimmed_{site}_{depth}_R1.fastq",
        R2 = "trimmed/trimmed_{site}_{depth}_R2.fastq"
    conda:
        "/faststorage/project/Sediment_Evolution/scripts/main_conda.yaml"
    shell: """
        mkdir -p trimmed
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2}
    """

rule gene_annotations:
    input:
        ref = "refs/{mag_origin}.{mag}.filtered.fa"
    output:
        out = "annotations/{mag_origin}.{mag}_eggnog.emapper.annotations",
        fna = "annotations/{mag_origin}.{mag}_eggnog.emapper.genepred.fasta"
    resources:
            cpus_per_task= 24,
            mem_mb= 32000
    conda:
        "/faststorage/project/Sediment_Evolution/scripts/eggnog_conda.yaml"
    shell: """
        emapper.py --cpu {resources.cpus_per_task} --itype genome -i {input.ref} --data_dir /home/eawi/eggnog-db/ --output annotations/$(basename {input.ref} .filtered.fa)_eggnog
    """

rule progigal:
    input:
        comb_fa = "combined_reference_{mag_origin}.fa",
    output:
        comb_genes = "combined_reference_{mag_origin}.fna"
    conda:
        "/faststorage/project/Sediment_Evolution/scripts/main_conda.yaml"
    shell: """
        prodigal -d {output.comb_genes} -i {input.comb_fa}
    """

rule bwa_index_mags:
    input:
        comb_fa = "combined_reference_{mag_origin}.fa",
    output:
        comb_fa_index = "combined_reference_{mag_origin}.fa.pac",
    shell:"""
    bwa index {input.comb_fa}
    """
rule first_read_mapping:
    input:
        comb_fa = "combined_reference_{mag_origin}.fa",
        comb_fa_index = "combined_reference_{mag_origin}.fa.pac",
        R1 = "trimmed/trimmed_" + my_site + "_{depth}_R1.fastq",
        R2 = "trimmed/trimmed_" + my_site + "_{depth}_R2.fastq"
    output:
        sam = "bams/temp_{mag_origin}_{depth}.sam"
    resources:
        cpus_per_task = 16,
        mem_mb = 16000
    conda:
        "/faststorage/project/Sediment_Evolution/scripts/main_conda.yaml"
    shell: """
        bwa mem -t {resources.cpus_per_task} {input.comb_fa} {input.R1} {input.R2} > {output.sam}
    """
rule sort_sam_file:
    input:
        sam = "bams/temp_{mag_origin}_{depth}.sam"
    output:
        bam = "bams/raw_{mag_origin}_{depth}.bam"
    resources:
        mem_mb = 8000,
        cpus_per_task = 4,
    shell: """
        samtools sort -o {output.bam} -O BAM -@ {resources.cpus_per_task} {input.sam} &&
        rm {input.sam}
    """

rule first_instrain:
    input:
        comb_fa = "combined_reference_{mag_origin}.fa",
        bam = "bams/raw_{mag_origin}_{depth}.bam",
        stb = "combined_reference_{mag_origin}.stb",
        comb_genes = "combined_reference_{mag_origin}.fna"
    params:
        isdir = "IS_raw_{mag_origin}_{depth}"  
    output:
        output = "IS_raw_{mag_origin}_{depth}/output/IS_raw_{mag_origin}_{depth}_genome_info.tsv",
        rdic = "IS_raw_{mag_origin}_{depth}/raw_data/Rdic.json",
    conda:
        "/faststorage/project/Sediment_Evolution/scripts/main_conda.yaml"
    resources:
        cpus_per_task = 32,
        mem_mb = 128000
    shell: """
       inStrain profile -p {resources.cpus_per_task} -s {input.stb} -g {input.comb_genes} -o {params.isdir} {input.bam} {input.comb_fa}
    """
#
rule filter_qnames:
    input:
        Rdic = "IS_raw_{mag_origin}_{depth}/raw_data/Rdic.json",
        stb = "combined_reference_{mag_origin}.stb",
        bam = "bams/raw_{mag_origin}_{depth}.bam",
    output:
        filtered = "bams/filtered_{mag_origin}_{depth}.bam",
    run:
        qnames = f"{my_site}_{wildcards.depth}_qnames.txt"
        with open(input.Rdic) as f:
            my_json = json.load(f)

        with open(input.stb) as f:
            contigs = [line.strip().split("\t")[0] for line in f]

        with open(qnames,"w") as f:
            for ctg in contigs:
                try:
                    for qname in set(q.split(" ")[0] for q in my_json[ctg].keys()):
                        f.write(qname)
                        f.write("\n")
                except KeyError as e:
                    print(ctg)
        os.system(f"samtools view -N {qnames} -O BAM -o {output.filtered} {input.bam}")

rule find_regions:
    input: 
        comb_fna = "combined_reference_{mag_origin}.fna"
    output:
        regions = "regions_{mag_origin}.tsv"
    run:
        region_list = []
        with open(input.comb_fna) as f:
            for line in f:
                if line[0] == ">":
                    ctg, start, end, direction, extra = line[1:-1].replace(" ","").split("#")
                    ctg = "_".join(ctg.strip().split("_")[:-1])
                    reg = f"{ctg}:{start}-{end}"
                    region_list.append(reg)
        with open(output.regions, "w") as f:
            for reg in region_list:
                f.write("\t".join([reg for i in range(len(my_depths))])+"\n") # repeat the region as many times as there are bams

rule normalize:
    input:
        filtered_bams = expand(["bams/filtered_{mo}_{d}.bam"], mo=mag_origin, d=my_depths),
        regions = f"regions_{mag_origin}.tsv"
    output:
        expand(["normalized_filtered_{mo}_{d}.bam"], mo=mag_origin, d=my_depths)
    
    shell: """
        /faststorage/project/Sediment_Evolution/scripts/normalize_paired_regions/normalize_paired_regions -r {input.regions} --bams {input.filtered_bams} && touch normalize_done.txt
    """

rule sort_and_index_normalized:
    input:
        normalized_bam = "bams/normalized_filtered_{mag_origin}_{depth}.bam"
    output:
        proc_bam = "bams/processed_{mag_origin}_{depth}.bam"
    conda:
        "/faststorage/project/Sediment_Evolution/scripts/main_conda.yaml"
    shell: """
        samtools sort -o {output.proc_bam} -O BAM {input.normalized_bam}
    """

rule second_pass_instrain:
    input:
        proc_bam = "bams/processed_{mag_origin}_{depth}.bam",
        stb = "combined_reference_{mag_origin}.stb",
        comb_genes = "combined_reference_{mag_origin}.fna",
        comb_fa = "combined_reference_{mag_origin}.fa"
    output:
        is_out_dir = directory("IS_proc_{mag_origin}_{depth}"),
        output = "IS_proc_{mag_origin}_{depth}/output/IS_proc_{mag_origin}_{depth}_genome_info.tsv"
    resources:
        cpus_per_task = 48,
        mem_mb = 64000
    conda:
        "/faststorage/project/Sediment_Evolution/scripts/main_conda.yaml"
    shell: """
       inStrain profile --skip_plot_generation -p {resources.cpus_per_task} -s {input.stb} -g {input.comb_genes} -o {output.is_out_dir} {input.proc_bam} {input.comb_fa}
    """
