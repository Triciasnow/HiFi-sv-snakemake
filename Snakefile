##
## Run with, for example:
# > snakemake --configfile config.yaml --cores 64 --use-conda --use-envmodules read_map_all
# Config can also be overwritten in the command line, for example:
# > snakemake --configfile config.yaml --config samples="Tomatoes" --use-conda --use-envmodules --cores 64 read_map_all
##

##
## Config:
##
configfile: "config.yaml"

REF=config['reference']
SAMPLES=config['samples'].split()
DEPTH=config['depth']
ALIGNERS=config['aligners'].split()

##
## main rules:
##

## mapped by mutiple aligner
rule read_map_all:
    input:
        expand('{sample}-{aligner}.sort.bam', sample=SAMPLES, aligner=['NGMLR','pbmm2']),
        expand('{sample}-{aligner}.{ext}', sample=SAMPLES, aligner=['minimap2'], ext=['bam']),
        expand('{sample}-{aligner}.{ext}', sample=SAMPLES, aligner=['nucmer'], ext=['delta'])

# calling by all software
rule call_original_all:
    input:
        expand('{sample}/{caller}/{sample}-{aligner}.{caller}.vcf', sample=SAMPLES, aligner=['NGMLR','pbmm2'], caller=['pbsv','SVIM','Sniffles','cuteSV']),
        expand('{sample}/{caller}/{sample}.{caller}.{ext}', sample=SAMPLES, caller=['Assemblytics'], ext=['bed']),
        expand('{sample}/{caller}/{sample}.{caller}.{ext}', sample=SAMPLES, caller=['SyRI'], ext=['vcf'])

rule call_filter_all:
    input:
        expand('{sample}/{type}/{sample}-{aligner}.{caller}.filter.vcf', sample=SAMPLES, aligner=['NGMLR','pbmm2'], type=['INS'], caller=['SVIM']),
        expand('{sample}/{type}/{sample}-{aligner}.{caller}.filter.vcf', sample=SAMPLES, aligner=['NGMLR','pbmm2'], type=['DEL'], caller=['Sniffles']),
        expand('{sample}/{type}/{sample}.{caller}.filter.vcf', sample=SAMPLES, type=['DEL'], caller=['Assemblytics']),
        expand('{sample}/{type}/{sample}-{aligner}.{caller}.filter.vcf', sample=SAMPLES, aligner=['NGMLR'], type=['INV'], caller=['Sniffles','cuteSV']),
        expand('{sample}/{type}/{sample}-{aligner}.{caller}.filter.vcf', sample=SAMPLES, aligner=['pbmm2'], type=['DUP'], caller=['pbsv','cuteSV']),
        expand('{sample}/{type}/{sample}-{aligner}.{caller}.filter.vcf', sample=SAMPLES, aligner=['NGMLR','pbmm2'], type=['TRA'], caller=['pbsv','Sniffles','cuteSV']),
        expand('{sample}/{type}/{sample}.{caller}.filter.vcf', sample=SAMPLES, type=['TRA'], caller=['SyRI'])

rule final_results_all:
    input:
        expand('{sample}/results/{type}.vcf', sample=SAMPLES, type=['INS','DEL','DUP','INV','TRA'])




##
## SVs calling
##

## Mapping samples to reference genome with NGMLR and pbmm2
# alignment with NGMLR
rule NGMLR_map:
    input:
        ref=expand('{genome}.fa',genome=REF),
        fq=expand('{sample}.fq.gz',sample=SAMPLES)
    output:
        '{sample}-NGMLR.bam'
    threads: 12
    log:
        'logs/{sample}-NGMLR.log.txt'
    conda: "./evn/ngmlr.yaml"
    params:
        para=config['reference']
    shell:
        '(ngmlr -r {input.ref} -q {input.fq} -t {threads} -x pacbio --rg-sm {wildcards.sample} --rg-lb {params.para} --rg-pl PacBio -o {output}) 2> {log}'
# alignment with pbmm2
rule pbmm2_map:
    input:
        ref=expand('{genome}.fa',genome=REF),
        fq=expand('{sample}.fq.gz',sample=SAMPLES)
    output:
        '{sample}-pbmm2.bam'
    params:
        para=config['reference']
    conda: "./evn/pbmm2.yaml"
    threads: 1
    log:
        'logs/{sample}-pbmm2.log.txt'
    shell:
        "(pbmm2 align {input} {output} --unmapped --preset CCS -N 2 --sort --rg '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{params.para}\\tPL:PB')' 2> {log}"


## index and sort mapped result
rule samtools_sort_index:
    input:
        bam_ngmlr=expand('{sample}-NGMLR.bam',sample=SAMPLES),
        bam_pbmm2=expand('{sample}-pbmm2.bam',sample=SAMPLES)
    output:
        sort_ngmlr='{sample}-NGMLR.sort.bam',
        sort_pbmm2='{sample}-pbmm2.sort.bam',
        index_ngmlr='{sample}-NGMLR.bam.bai',
        index_pbmm2='{sample}-pbmm2.bam.bai'
    conda: "./evn/samtools.yaml"
    threads: 4
    log:
        'logs/{sample}-sort_index.log.txt'
    shell:
        '(samtools sort {input.bam_ngmlr} -o {output.sort_ngmlr} | samtools index) 2> {log}' \
        '(samtools sort {input.bam_pbmm2} -o {output.sort_pbmm2} | samtools index) 2>> {log}'



## Mapping assembly genome to reference genome with MuMmer4,minimap2,and LASTZ
# alignment with MuMmer4:
rule MuMmer4_map:
    input:
        ref=expand('{genome}.fa',genome=REF),
        asm=expand('{sample}.fa',sample=SAMPLES)
    output:
        '{sample}-nucmer.delta'
    threads: 12
    log:
        'logs/{sample}-nucmer.log.txt'
    params:
        para=config['nucmer_para'],
        name=config['reference']
    conda: "./evn/mummer4.yaml"
    shell:
        '(nucmer --maxmatch {params.para} {input.ref} {input.asm} -p {wildcards.sample}-{params.name} -t {threads}) 2> {log}' \
        '(mv {wildcards.sample}-{params.name}.delta {output})'

# alignment with minimap2:
rule minimap2_map:
    input:
        ref=expand('{genome}.fa',genome=REF),
        asm=expand('{sample}.fa',sample=SAMPLES)
    output:
        sam=temp('{sample}-minimap2.sam'),
        bam='{sample}-minimap2.bam'
    conda: "./evn/minimap2.yaml"
    threads: 2
    log:
        'logs/{sample}-minimap2.log.txt'
    shell:
        '(minimap2 -ax asm5 --eqx {input.ref} {input.asm} > {output.sam} | samtools view -b -o {output.bam}) >2 {log}'


## sv calling by software which based on reads alignment
# calling by pbsv
rule pbsv_call:
    input:
        ref=expand('{genome}.fa',genome=REF),
        bam=expand('{sample}-{aligner}.sort.bam',sample=SAMPLES,aligner=['NGMLR','pbmm2']),
        index='{sample}-{aligner}.bam.bai'
    output:
        repeat='{sample}/pbsv/{sample}-{aligner}.trf.bed',
        vcf='{sample}/pbsv/{sample}-{aligner}.pbsv.vcf',
        svsig = temp("{sample}/pbsv/{sample}-{aligner}.pbsv.svsig.gz"),
    threads: 2
    conda: "./evn/pbsv.yaml"
    params:
        para=config['trf_para'],
        depth=config['depth']
    log:
        tdm='logs/{sample}-{aligner}-tandem.log.txt',
        sv='logs/{sample}-{aligner}.pbsv.log.txt',
    shell:
        '(trf {input.ref} {params.para}) 2> {log.tdm}' \
        '(pbsv discover {input.bam} --tandem-repeats {output.repeat} {output.svsig}) 2> {log.sv}' \
        '(pbsv call --ccs -x {params.depth} {input.ref} {output.svsig} {output.vcf}) 2>> {log.sv}'

# calling by SVIM
rule SVIM_call:
    input:
        ref=expand('{genome}.fa',genome=REF),
        bam=expand('{sample}-{aligner}.sort.bam',sample=SAMPLES,aligner=['NGMLR','pbmm2'])
    output:
        '{sample}/SVIM/{sample}-{aligner}.SVIM.vcf'
    threads: 2
    log:
        'logs/{sample}-{aligner}.SVIM.log.txt'
    params:
        outdir='{sample}/SVIM'
    conda: "./evn/svim.yaml"
    shell:
        '(svim alignment --sequence_alleles --read_names {input.ref} {input.bam} {params.outdir}/) 2> {log}' \
        '(mv {params.outdir}/final_results.vcf {output})'

# calling by Sniffles
if config['depth'] >= 30:
    rule Sniffles_call:
        input:
            bam=expand('{sample}-{aligner}.sort.bam',sample=SAMPLES,aligner=['NGMLR','pbmm2'])
        output:
            '{sample}/Sniffles/{sample}-{aligner}.Sniffles.vcf'
        threads: 2
        log:
            'logs/{sample}-{aligner}.Sniffles.log.txt'
        params:
            para=int(config['depth'])/10
        conda: "./evn/sniffles.yaml"
        shell:
            '(sniffles --skip_parameter_estimation -s {params.para} --ccs_reads -m {input} -v {output}) 2> {log}'
else:
    rule Sniffles_call:
        input:
            bam=expand('{sample}-{aligner}.sort.bam',sample=SAMPLES,aligner=['NGMLR','pbmm2'])
        output:
            '{sample}/Sniffles/{sample}-{aligner}.Sniffles.vcf'
        threads: 2
        log:
            'logs/{sample}-{aligner}.Sniffles.log.txt'
        params:
            para=2
        conda: "./evn/sniffles.yaml"
        shell:
            '(sniffles --skip_parameter_estimation -s {params.para} --ccs_reads -m {input} -v {output}) 2> {log}'
# calling by cuteSV
rule cuteSV_call:
    input:
        ref=expand('{genome}.fa',genome=REF),
        bam=expand('{sample}-{aligner}.sort.bam',sample=SAMPLES,aligner=['NGMLR','pbmm2'])
    output:
        '{sample}/cuteSV/{sample}-{aligner}.cuteSV.vcf'
    threads: 2
    log:
        'logs/{sample}-{aligner}.cuteSV.log.txt'
    params:
        map_name=config['aligners']
        work_dir='{sample}/cuteSV/{sample}-{aligner}_tmp'
    conda: "./evn/cutesv.yaml"
    shell:
        'mkdir {params.work_dir}' \
        '(cuteSV {input} {output} {params.work_dir} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --threads {threads} --sample {wildcards.sample}-{params.map_name} --min_support 10 --genotype) 2> {log}'

## sv calling by software which based on genoem alignment
# calling by Assemblytics
rule Assemblytics_call:
    input:
        '{sample}-nucmer.delta'
    output:
        '{sample}/Assemblytics/{sample}.Assemblytics.bed'
    threads: 2
    log:
        'logs/{sample}.Assemblytics.log.txt'
    params:
        name=config['reference']
    conda: "./evn/assemblytics.yaml"
    shell:
        '(gzip {input} | scripts/Assemblytics {wildcards.sample}-{params.name} 10000 20 100000) 2> {log}' \
        'mv {wildcards.sample}-{params.name}/{wildcards.sample}-{params.name}.Assemblytics_structural_variants.bed {output}'
# calling by SyRI
rule SyRI_call:
    input:
        ref=expand('{genome}.fa',genome=REF),
        bam='{sample}-minimap2.bam',
        asm=expand('{sample}.fa',sample=SAMPLES)
    output:
        ref_chr='{genome}.rechr.fa',
        asm_chr='{sample}.rechr.fa',
        vcf='{sample}/SyRI/{sample}.SyRI.vcf'
    threads: 2
    conda: "./evn/syri.yaml"
    params:
        max_chr=config['Chrs_ref']
    log:
        'logs/{sample}.SyRI.log.txt'
    shell:
        'python ref_chr.py -i {input.ref} -o {output.ref_chr} -c {params.max_chr}' \
        'python asm_chr.py -i {input.asm} -o {output.asm_chr}' \
        '(syri -c {input.bam} -r {output.ref_chr} -q {output.asm_chr} -k -F B --nosnp) 2> {log}' \
        'mv syri.vcf {output.vcf}'


##
## Final SVs set generate
##


## INS detection: pbmm2+NGMLR-SVIM
rule INS_call_filter:
    input:
        '{sample}/SVIM/{sample}-{aligner}.SVIM.vcf'
    output:
        '{sample}/INS/{sample}-{aligner}.SVIM.filter.vcf'
    threads: 1
    params:
        qual=config['svim_qual']
    shell:
        'python svim_ins_filter.py -i {input} -q {params.qual} -o {output}'
# INS merge:
rule INS_merge:
    input:
        '{sample}/INS/{sample}-pbmm2.SVIM.filter.vcf',
        '{sample}/INS/{sample}-NGMLR.SVIM.filter.vcf'
    output:
        '{sample}/results/INS.vcf'
    log:
        'logs/{sample}.INS.log.txt'
    threads: 1
    conda: "./evn/survivor.yaml"
    shell:
        '(ls {input} > sample_files_ins)' \
        '(SURVIVOR merge sample_files_ins 1000 1 1 1 0 20 {output}) 2> {log}'

##
# DEL detection: pbmm2+NGMLR-Sniffles; Assemblytics
rule DEL_Assemblytics_format:
    input:
        bed='{sample}/Assemblytics/{sample}.Assemblytics.bed',
        ref=expand('{genome}.fa',genome=REF),
        txt=expand('{anno}.txt',anno=config['assemblytics_anno'])
    output:
        vcf_temp=temp('{sample}/DEL/{sample}.Assemblytics.temp.vcf'),
        norm=temp('{sample}/DEL/{sample}.Assemblytics.norm.vcf'),
        vcf='{sample}/DEL/{sample}.Assemblytics.filter.vcf'
    threads: 1
    shell:
        'python assemblytics_format_change.py -i {input.bed} -a {input.txt} -o {output.vcf_temp}' \
        'bcftools norm -f {input.ref} -m-both -c s -o {output.norm} {output.vcf_temp}' \
        'python genotype.py -i {output.vcf_temp} -o {output.vcf}'
rule DEL_Sniffles_filter:
    input:
        '{sample}/Sniffles/{sample}-{aligner}.Sniffles.vcf'
    output:
        '{sample}/DEL/{sample}-{aligner}.Sniffles.filter.vcf'
    threads: 1
    run:
        shell('python sniffles_del_fileter.py -i {input} -o {output}')
rule DEL_merge:
    input:
        sniffles_ngmlr='{sample}/DEL/{sample}-NGMLR.Sniffles.filter.vcf',
        sniffles_pbmm2='{sample}/DEL/{sample}-pbmm2.Sniffles.filter.vcf',
        Assemblytics='{sample}/DEL/{sample}.Assemblytics.filter.vcf'
    output:
        '{sample}/results/DEL.vcf'
    threads: 1
    log:
        'logs/{sample}.DEL.log.txt'
    run:
        shell('(ls {input} > sample_files_del)')
        shell('(SURVIVOR merge sample_files_del 1000 1 1 1 0 20 {output}) 2> {log}')


##
# INV detection: NGMLR-cuteSV,Sniffles
rule INV_call_filter:
    input:
        cuteSV='{sample}/cuteSV/{sample}-NGMLR.cuteSV.vcf',
        Sniffles='{sample}/Sniffles/{sample}-NGMLR.Sniffles.vcf'
    output:
        cuteSV_vcf='{sample}/INV/{sample}-NGMLR.cuteSV.filter.vcf',
        Sniffles_vcf='{sample}/INV/{sample}-NGMLR.Sniffles.filter.vcf'
    threads: 1
    run:
        shell('python sniffles_inv_fileter.py -i {input.Sniffles} -o {output.Sniffles_vcf}')
        shell('python cutesv_inv_fliter.py -i {input.cuteSV} -o {output.cuteSV_vcf}')
rule INV_merge:
    input:
        cuteSV_vcf='{sample}/INV/{sample}-NGMLR.cuteSV.filter.vcf',
        Sniffles_vcf='{sample}/INV/{sample}-NGMLR.Sniffles.filter.vcf'
    output:
        '{sample}/results/INV.vcf'
    threads: 1
    log:
        'logs/{sample}.INV.log.txt'
    run:
        shell('(ls {input} > sample_files_inv)')
        shell('(SURVIVOR merge sample_files_inv 1000 1 1 1 0 20 {output}) 2> {log}')


# DUP detection: pbmm2-pbsv,cuteSV
rule DUP_call_filter:
    input:
        pbsv='{sample}/pbsv/{sample}-pbmm2.pbsv.vcf',
        cuteSV='{sample}/cuteSV/{sample}-pbmm2.cuteSV.vcf'
    output:
        pbsv_vcf='{sample}/DUP/{sample}-pbmm2.pbsv.filter.vcf',
        cuteSV_vcf='{sample}/DUP/{sample}-pbmm2.cuteSV.filter.vcf'
    threads: 1
    run:
        shell('python pbsv_dup_fileter.py -i {input.pbsv} -o {output.pbsv_vcf}')
        shell('python cutesv_dup_fliter.py -i {input.cuteSV} -o {output.cuteSV_vcf}')
rule DUP_merge:
    input:
        pbsv_vcf='{sample}/DUP/{sample}-pbmm2.pbsv.filter.vcf',
        cuteSV_vcf='{sample}/DUP/{sample}-pbmm2.cuteSV.filter.vcf'
    output:
        '{sample}/results/DUP.vcf'
    threads: 1
    log:
        'logs/{sample}.DUP.log.txt'
    run:
        shell('(ls {input} > sample_files_dup)')
        shell('(SURVIVOR merge sample_files_dup 1000 1 1 1 0 20 {output}) 2> {log}')


# TRA detection: pbmm2+NGMLR-pbsv,Sniffles,cuteSV; SyRI
rule TRA_SyRI_format:
    input:
        '{sample}/SyRI/{sample}.SyRI.vcf'
    output:
        '{sample}/TRA/{sample}.SyRI.filter.vcf'
    threads: 1
    run:
        shell('python syri_tra_filter.py -i {input} -o {output}')
rule TRA_filter:
    input:
        pbsv='{sample}/pbsv/{sample}-{aligner}.pbsv.vcf',
        Sniffles='{sample}/Sniffles/{sample}-{aligner}.Sniffles.vcf',
        cuteSV='{sample}/cuteSV/{sample}-{aligner}.cuteSV.vcf'
    output:
        pbsv_vcf='{sample}/TRA/{sample}-{aligner}.pbsv.filter.vcf',
        Sniffles_vcf='{sample}/TRA/{sample}-{aligner}.Sniffles.filter.vcf',
        cuteSV_vcf='{sample}/TRA/{sample}-{aligner}.cuteSV.filter.vcf'
    threads: 1
    run:
        shell('python pbsv_tra_filter.py -i {input.pbsv} -o {output.pbsv_vcf}')
        shell('python sniffles_tra_filter.py -i {input.Sniffles} -o {output.Sniffles_vcf}')
        shell('python cutesv_tra_filter.py -i {input.cuteSV} -o {output.cuteSV_vcf}')
rule TRA_merge:
    input:
        pbsv_vcf1='{sample}/TRA/{sample}-pbmm2.pbsv.filter.vcf',
        pbsv_vcf2='{sample}/TRA/{sample}-NGMLR.pbsv.filter.vcf',
        Sniffles_vcf1='{sample}/TRA/{sample}-pbmm2.Sniffles.filter.vcf',
        Sniffles_vcf2='{sample}/TRA/{sample}-NGMLR.Sniffles.filter.vcf',
        cuteSV_vcf1='{sample}/TRA/{sample}-pbmm2.cuteSV.filter.vcf',
        cuteSV_vcf2='{sample}/TRA/{sample}-NGMLR.cuteSV.filter.vcf',
        SyRI_vcf='{sample}/TRA/{sample}.SyRI.filter.vcf'
    output:
        '{sample}/TRA/TRA_temp.vcf'
    threads: 1
    log:
        'logs/{sample}.TRA.log.txt'
    run:
        shell('(ls {input} > sample_files_tra)')
        shell('(SURVIVOR merge sample_files_tra 1000 2 1 1 0 20 {output}) 2> {log}')
rule TRA_remove_repeat:
    input:
        '{sample}/TRA/TRA_temp.vcf'
    output:
        '{sample}/results/TRA.vcf'
    threads: 1
    run:
        shell('python tra_remove_repeat.py -i {input} -o {output}')
