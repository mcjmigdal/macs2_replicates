# input
SAMPLES = ["DSHBco.filtered", "DSHB39.filtered", "Isl.filtered"]
CONTROL = ["input.filtered.bam"]

# params
cutoff = 2
minlen = 600
maxgap = 100

rule call_peaks_get_bdg:
    input:
      treatment = "{sample}.bam",
      control = expand("{control}", control=CONTROL)
    output:
      "{sample}"
    shell:
      "macs2 callpeak --treatment {input.treatment} --control {input.control} "
      "--name {output} --outdir . --format BAM --gsize 1.37e9 --qvalue 0.01 --bdg "
      "--nomodel --extsize 600 2> {output}.log"

rule get_pvalues:
    input:
      pileup = "{sample}_treat_pileup.bdg",
      localbias = "{sample}_control_lambda.bdg"
    output:
      "{sample}.pvalue.bdg"
    shell:
      "macs2 bdgcmp -t {input.pileup} -c {input.localbias} -m ppois -o {output}"

rule combine_pvalues:
    input:
        expand("{sample}.pvalue.bdg", sample=SAMPLES)
    output:
        "combined.qvalue.bdg"
    shell:
        "Rscript combine_pvalues.R {output} {input} 2> /dev/null"

rule call_peaks:
    input:
        "combined.qvalue.bdg"
    output:
        "combined.narrowPeak"
    params:
        cutoff = cutoff, # -log10(1e-5)
        minlen = int(minlen), # minimum length of peak
        maxgap = int(maxgap) # maximum gap between significant points *in a peak*
    shell:
        "macs2 bdgpeakcall -i {input} -c {params.cutoff} -l {params.minlen} -g {params.maxgap} -o {output}"
        
# rule clean working directory?
