# input
SAMPLES = ["kdrl", "fabp", "hand"]
REPLICATES = ["1", "2", "3"]

# params
max_fragment_length = 128
gsize = 1.37e9
keepduplicates = "auto"
atac_tag_extsize = 50 # both sides
slocal_extsize = 2500
llocal_extsize = 5000

# peak detection params
cutoff = 5 # -log10(1e-5)
minlen = 100 # minimum length of peak
maxgap = 75 # maximum gap between significant points *in a peak*

rule filter_shift_bam:
    input:
        "{sample}_{replicate}.bam"
    output:
        "{sample}_{replicate}.shifted.bam"
    params:
        th = max_fragment_length
    shell:
      "alignmentSieve -b {input} -o {output} --ATACshift --maxFragmentLength {params.th}"
  

rule filter_duplicates:
    input:
        "{sample}_{replicate}.shifted.bam"
    output:
        out = "{sample}_{replicate}.filterdup.bed",
        log = "{sample}_{replicate}.log"
    params:
        format = "BAM",
        gsize = int(gsize),
        keepduplicates = keepduplicates,
        tsize = 1
    shell:
        "macs2 filterdup -i {input} -f {params.format} --gsize {params.gsize} "
        "--tsize {params.tsize} --keep-dup {params.keepduplicates} -o {output.out} 2>> {output.log}"
        
rule generate_coverage_track:
    input:
        "{sample}_{replicate}.filterdup.bed"
    output:
        "{sample}_{replicate}.filterdup.pileup.bdg"
    params:
        format = "BED",
        extsize = int(atac_tag_extsize)
    shell:
        "macs2 pileup -i {input} -B -f {params.format} --extsize {params.extsize} -o {output} 2> /dev/null"
        
rule build_slocal_bias_track:
    input:
        "{sample}_{replicate}.filterdup.bed"
    output:
        slocal = "{sample}_{replicate}.slocal.bdg",
        slocalnorm = "{sample}_{replicate}.slocal.norm.bdg"
    params:
        format = "BED",
        extsize = int(slocal_extsize),
        scale = atac_tag_extsize / float(slocal_extsize) # d / slocal normalize for the coverage track
    shell:
        "macs2 pileup -i {input} -B -f {params.format} --extsize {params.extsize} -o {output.slocal} 2> /dev/null;"
        "macs2 bdgopt -i {output.slocal} -m multiply -p {params.scale} -o {output.slocalnorm} 2> /dev/null"

rule build_llocal_bias_track:
    input:
        "{sample}_{replicate}.filterdup.bed"
    output:
        llocal = "{sample}_{replicate}.llocal.bdg",
        llocalnorm = "{sample}_{replicate}.llocal.norm.bdg"
    params:
        format = "BED",
        extsize = int(llocal_extsize),
        scale = atac_tag_extsize / float(llocal_extsize) # d / llocal normalize for the coverage track
    shell:
        "macs2 pileup -i {input} -B -f {params.format} --extsize {params.extsize} -o {output.llocal} 2> /dev/null;"
        "macs2 bdgopt -i {output.llocal} -m multiply -p {params.scale} -o {output.llocalnorm} 2> /dev/null"

rule generate_maximum_background_noise:
    input:
        slocal = "{sample}_{replicate}.slocal.norm.bdg",
        llocal = "{sample}_{replicate}.llocal.norm.bdg",
        log = "{sample}_{replicate}.log"
    output:
        normlocal = "{sample}_{replicate}.local.raw.tmp.bdg",
        localbias = "{sample}_{replicate}.local.raw.bdg"
    params:
        genomebck = "%.10f" % ((2 * atac_tag_extsize) / gsize) # the_number_of_control_reads*fragment_length/genome_size -- this might be substituded with pseudocounts as it might fail if the result genomebck is 0
    shell:
        "genomebck=$(grep 'tags after filtering in alignment file' {input.log} | rev | cut -d : -f 1 | rev);"
        "genomebck=$(echo $genomebck \* {params.genomebck} | bc -l);"
        "echo genome wide backgroud: $genomebck >> {input.log};"
        "macs2 bdgcmp -m max -t {input.slocal} -c {input.llocal} -o {output.normlocal} 2> /dev/null;"
        "macs2 bdgopt -i {output.normlocal} -m max -p $genomebck -o {output.localbias} 2> /dev/null"

rule get_pvalues:
    input:
        pileup = "{sample}_{replicate}.filterdup.pileup.bdg",
        localbias = "{sample}_{replicate}.local.raw.bdg"
    output:
        "{sample}_{replicate}.pvalue.bdg"
    shell:
        "macs2 bdgcmp -t {input.pileup} -c {input.localbias} -m ppois -o {output}"
        
rule combine_pvalues:
    input:
        expand("{sample}_{replicate}.pvalue.bdg", sample=SAMPLES, replicate=REPLICATES)
    output:
        "{sample}.combined.qvalue.bdg"
    shell:
        "Rscript combine_pvalues.R {output} {input} 2> /dev/null"

rule call_peaks:
    input:
        "{sample}.combined.qvalue.bdg"
    output:
        "{sample}.narrowPeak"
    params:
        cutoff = cutoff, # -log10(1e-5)
        minlen = int(minlen), # minimum length of peak
        maxgap = int(maxgap) # maximum gap between significant points *in a peak*
    shell:
        "macs2 bdgpeakcall -i {input} -c {params.cutoff} -l {params.minlen} -g {params.maxgap} -o {output}"
        
# rule clean working directory?
