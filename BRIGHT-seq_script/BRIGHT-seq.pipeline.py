#!/usr/bin/env runsnakemake

################################################################################
# Pipeline: BRIGHT-seq Analysis Workflow
# Description: Comprehensive pipeline for RNA editing and methylation analysis
#              Processes paired-end sequencing data through QC, alignment,
#              deduplication, and variant calling with strand-specific analysis
# Author: Ziang Lu
# Date: 2025.10.13
################################################################################

# Load configuration file
configfile: "BRIGHT-seq.pipeline.config.yaml"

################################################################################
# Configuration variables
################################################################################

samples = config["samples"]          # All sample names
ipsamples = config["ipsamples"]      # IP (immunoprecipitation) samples
insamples = config["insamples"]      # Input control samples
refs = config["refs"]                # Reference chromosomes/contigs
ori = config["ori"]                  # Strand orientation (F1-R2, R1-F2)
type = config["muttype"]             # Mutation type (c2t, g2a, etc.)
nowref = config["myref"]             # Current reference genome
datadir = config["datadir"]          # Input data directory
outdir = config["outdir"]            # Output directory
refdir = config["refdir"]            # Reference files directory

################################################################################
# Target rule - defines all final outputs
################################################################################

rule all:
    input:
        # Step 1: Quality control - trimmed reads
        expand(outdir + "/trimgalore/{sample}_R1_val_1.fq.gz", sample=samples),
        expand(outdir + "/trimgalore/{sample}_R2_val_2.fq.gz", sample=samples),
        
        # Step 2: Bowtie2 alignment to reference genome
        expand(outdir + "/bowtie2_{nowref}/{sample}.aligned.bam", sample=samples, nowref=nowref),
        
        # Step 3: BAM filtering and strand separation
        expand(outdir + "/bowtie2_{nowref}/{sample}.tmp2.bam", sample=samples, nowref=nowref),
        expand(outdir + "/bowtie2_{nowref}/{sample}.tmp2.bai", sample=samples, nowref=nowref),
        expand(outdir + "/bowtie2_{nowref}/{sample}.R1posR2neg.bam", sample=samples, nowref=nowref),
        expand(outdir + "/bowtie2_{nowref}/{sample}.R1negR2pos.bam", sample=samples, nowref=nowref),
        
        # Step 4: Coverage extraction
        expand(outdir + "/bowtie2_{nowref}/{sample}.aligned.cov", sample=samples, nowref=nowref),
        expand(outdir + "/bowtie2_{nowref}/{sample}.R1posR2neg.cov", sample=samples, nowref=nowref),
        expand(outdir + "/bowtie2_{nowref}/{sample}.R1negR2pos.cov", sample=samples, nowref=nowref),
        
        # Step 5: HISAT-3N bisulfite alignment
        expand(outdir + "/hisat3n/{sample}_{type}.sorted.F1-R2.bam", sample=samples, type=type),
        expand(outdir + "/hisat3n/{sample}_{type}.sorted.R1-F2.bam", sample=samples, type=type),
        
        # Step 6: Remove duplicates
        expand(outdir + "/hisat3n/{sample}_{type}.rmdup.{ori}.bam", sample=samples, ori=ori, type=type),
        
        # Step 7-8: REDItools2 variant calling and counting
        expand(outdir + "/reditools2/{sample}_{type}.redi.{ori}.{ref}.txt.gz", sample=samples, ori=ori, ref=refs, type=type),
        expand(outdir + "/reditools2/{sample}_{type}.redi.{ori}.{ref}.allmut.txt", sample=samples, ori=ori, ref=refs, type=type),
        expand(outdir + "/reditools2/{sample}_{type}.redi.{ori}.{ref}.fourbase.txt", sample=samples, ori=ori, ref=refs, type=type),
        
        # Step 9-14: SNP/mutation extraction and annotation
        expand(outdir + "/snp/{insample}_{type}.redi.{ori}.{ref}.snp.txt", insample=insamples, ori=ori, ref=refs, type=type),
        expand(outdir + "/mut/{ipsample}_{type}.redi.{ori}.{ref}.allmut.txt", ipsample=ipsamples, ori=ori, ref=refs, type=type),
        expand(outdir + "/snp/{insample}_{type}.redi.{ori}.newsnp.site.txt", insample=insamples, ori=ori, type=type),
        expand(outdir + "/snp/redi.{ori}.allsnp.site.txt", ori=ori),
        expand(outdir + "/mut/{ipsample}_{type}.redi.{ori}.mut.txt", ipsample=ipsamples, ori=ori, type=type),
        expand(outdir + "/mod/{ipsample}_{type}.redi.mod.anno.txt", ipsample=ipsamples, type=type),

################################################################################
# Step 1: Quality control with Trim Galore
################################################################################

rule my_step_01_trimgalore_ex:
    input:
        fq1 = datadir + "/{sample}_R1.fq.gz",
        fq2 = datadir + "/{sample}_R2.fq.gz",
    output:
        fq1 = outdir + "/trimgalore/{sample}_R1_val_1.fq.gz",
        fq2 = outdir + "/trimgalore/{sample}_R2_val_2.fq.gz",
    params:
        sedoutdir = outdir + "/trimgalore/"
    log:
        outdir + "/trimgalore/{sample}.trimgalore.mylog"
    threads: 4
    shell:
        """
        trimgalore.normal.sh {input.fq1} {input.fq2} {params.sedoutdir} {threads} &> {log}
        """

################################################################################
# Step 2: Bowtie2 alignment to reference genome
################################################################################

rule my_step_02_bowtie2_ref:
    input:
        fq1 = outdir + "/trimgalore/{sample}_R1_val_1.fq.gz",
        fq2 = outdir + "/trimgalore/{sample}_R2_val_2.fq.gz",
    output:
        outdir + "/bowtie2_{nowref}/{sample}.aligned.bam",
    params:
        genome = refdir + "/bowtie2_index_{nowref}/bowtie2_index_{nowref}",
    log:
        outdir + "/bowtie2_{nowref}/{sample}.mylog"
    threads: 4
    shell:
        """
        bowtie2_mismatch_N_1.sh {threads} {params.genome} {input.fq1} {input.fq2} {output} &> {log}
        """

################################################################################
# Step 3: BAM filtering and strand-specific splitting
################################################################################

rule my_step_03_bamextract:
    input:
        outdir + "/bowtie2_{nowref}/{sample}.aligned.bam",
    output:
        tmpbam2 = outdir + "/bowtie2_{nowref}/{sample}.tmp2.bam",
        tmpbai2 = outdir + "/bowtie2_{nowref}/{sample}.tmp2.bai",
        pos_strand = outdir + "/bowtie2_{nowref}/{sample}.R1posR2neg.bam",
        neg_strand = outdir + "/bowtie2_{nowref}/{sample}.R1negR2pos.bam",
    params:
        pos_bai = outdir + "/bowtie2_{nowref}/{sample}.R1posR2neg.bai",
        neg_bai = outdir + "/bowtie2_{nowref}/{sample}.R1negR2pos.bai",
        tmpbam1 = outdir + "/bowtie2_{nowref}/{sample}.tmp1.bam",
    log:
        outdir + "/bowtie2_{nowref}/{sample}.bamextract.mylog"
    threads: 1
    shell:
        """
        # Filter: primary alignment, proper pair, mapped, mapQ > 1
        bamtools filter -in {input} -isPrimaryAlignment true -isProperPair true \
            -isMapped true -mapQuality '>1' > {params.tmpbam1} 2> {log}
        
        # Sort and index filtered BAM
        samtools sort {params.tmpbam1} -o {output.tmpbam2} &>> {log}
        samtools index {output.tmpbam2} {output.tmpbai2}
        
        # Split by strand orientation
        python split_bam_specific_strand.py {output.tmpbam2} {output.pos_strand} {output.neg_strand}
        samtools index {output.pos_strand} {params.pos_bai}
        samtools index {output.neg_strand} {params.neg_bai}
        
        # Clean up temporary files
        rm -f {params.tmpbam1}
        """

################################################################################
# Step 4: Extract base coverage for all BAM files
################################################################################

rule my_step_04_cov_extract:
    input:
        align = outdir + "/bowtie2_{nowref}/{sample}.tmp2.bam",
        extra = outdir + "/bowtie2_{nowref}/{sample}.R1posR2neg.bam",
        extra2 = outdir + "/bowtie2_{nowref}/{sample}.R1negR2pos.bam",
    output:
        align_cov = outdir + "/bowtie2_{nowref}/{sample}.aligned.cov",
        extra_cov = outdir + "/bowtie2_{nowref}/{sample}.R1posR2neg.cov",
        extra_cov2 = outdir + "/bowtie2_{nowref}/{sample}.R1negR2pos.cov",
    params:
        fasta = refdir + "/fasta/{nowref}.fasta",
        tmpcov = outdir + "/bowtie2_{nowref}/{sample}.tmp.cov",
    log:
        outdir + "/bowtie2_{nowref}/{sample}.covextract.mylog"
    threads: 1
    shell:
        """
        # Extract coverage for main alignment
        python pysam_allbase_cov.py {input.align} > {params.tmpcov} 2>> {log}
        basecov_compare_ref.pl {params.fasta} {params.tmpcov} > {output.align_cov} 2>> {log}
        
        # Extract coverage for positive strand
        python pysam_allbase_cov.py {input.extra} > {params.tmpcov} 2>> {log}
        basecov_compare_ref.pl {params.fasta} {params.tmpcov} > {output.extra_cov} 2>> {log}
        
        # Extract coverage for negative strand
        python pysam_allbase_cov.py {input.extra2} > {params.tmpcov} 2>> {log}
        basecov_compare_ref.pl {params.fasta} {params.tmpcov} > {output.extra_cov2} 2>> {log}
        
        # Clean up temporary file
        rm -f {params.tmpcov}
        """

################################################################################
# Step 5: HISAT-3N alignment for bisulfite-treated reads
################################################################################

rule my_step_05_hisat3n:
    input:
        fq1 = outdir + "/trimgalore/{sample}_R1_val_1.fq.gz",
        fq2 = outdir + "/trimgalore/{sample}_R2_val_2.fq.gz",
    output:
        pos = outdir + "/hisat3n/{sample}_{type}.sorted.F1-R2.bam",
        neg = outdir + "/hisat3n/{sample}_{type}.sorted.R1-F2.bam",
    params:
        genome = refdir + "/hisat-3n_index_hg38.ucsc.analysisSet_{type}/hg38_{type}",
        sedoutdir = outdir + "/hisat3n/",
        sam = outdir + "/hisat3n/{sample}_{type}.sam",
        bam = outdir + "/hisat3n/{sample}_{type}.bam",
        bai = outdir + "/hisat3n/{sample}_{type}.sorted.bai",
        sortbam = outdir + "/hisat3n/{sample}_{type}.sorted.bam",
        pos_bai = outdir + "/hisat3n/{sample}_{type}.sorted.F1-R2.bai",
        neg_bai = outdir + "/hisat3n/{sample}_{type}.sorted.R1-F2.bai",
    log:
        outdir + "/hisat3n/{sample}_{type}.mylog",
    threads: 16
    shell:
        """
        cd {params.sedoutdir}
        
        # HISAT-3N alignment
        hisat3n.normal.sh {params.genome} {input.fq1} {input.fq2} {params.sam} {wildcards.type} {threads} 2> {log}
        
        # Convert SAM to BAM and sort
        samtools view -b -o {params.bam} -@ {threads} {params.sam} 2>> {log}
        samtools sort {params.bam} -o {params.sortbam} -@ {threads} 2>> {log}
        samtools index -@ {threads} {params.sortbam} {params.bai} 2>> {log}
        
        # Split by strand orientation
        python split_bam_specific_strand.py {params.sortbam} {output.pos} {output.neg} 2>> {log}
        samtools index -@ {threads} {output.pos} {params.pos_bai} 2>> {log}
        samtools index -@ {threads} {output.neg} {params.neg_bai} 2>> {log}
        
        # Clean up temporary files
        rm -f {params.bam} {params.sortbam}
        """

################################################################################
# Step 6: Remove PCR duplicates with Picard
################################################################################

rule my_step_06_rmdup:
    input:
        outdir + "/hisat3n/{sample}_{type}.sorted.{ori}.bam",
    output:
        outdir + "/hisat3n/{sample}_{type}.rmdup.{ori}.bam",
    params:
        tmpbam1 = outdir + "/hisat3n/{sample}_{type}.tmp1.{ori}.bam",
        rmdup = outdir + "/hisat3n/{sample}_{type}.rmdup.{ori}.txt",
        bai = outdir + "/hisat3n/{sample}_{type}.rmdup.{ori}.bai",
        sam = outdir + "/hisat3n/{sample}_{type}.sam",
    log:
        outdir + "/hisat3n/{sample}.rmdup.mylog"
    threads: 8
    shell:
        """
        # Remove duplicates with Picard
        removedup.picard.sh {input} {params.tmpbam1} {params.rmdup} &> {log}
        
        # Sort and index deduplicated BAM
        samtools sort -o {output} {params.tmpbam1} &>> {log}
        rm -f {params.bai}
        samtools index {output} {params.bai}
        
        # Clean up temporary files
        rm -f {params.tmpbam1} {params.sam}
        """

################################################################################
# Step 7: REDItools2 variant calling
################################################################################

rule my_step_07_REDItools2:
    input:
        outdir + "/hisat3n/{sample}_{type}.rmdup.{ori}.bam",
    output:
        outdir + "/reditools2/{sample}_{type}.redi.{ori}.{ref}.txt.gz",
    params:
        fasta = refdir + "/fasta/hg38.ucsc.analysisSet.fa",
        out = outdir + "/reditools2/{sample}_{type}.redi.{ori}.{ref}.txt",
    log:
        outdir + "/reditools2/{sample}.table.log"
    threads: 1
    shell:
        """
        # Run REDItools2 for variant calling
        python /home/reditools2.0/src/cineca/reditools.py \
            -f {input} -o {params.out} -r {params.fasta} -g {wildcards.ref} 2> {log}
        
        # Compress output
        gzip -c {params.out} 1> {output} 2>> {log}
        rm -f {params.out}
        """

################################################################################
# Step 8: Extract mutation and base counts from REDItools2 output
################################################################################

rule my_step_08_count:
    input:
        outdir + "/reditools2/{sample}_{type}.redi.{ori}.{ref}.txt.gz",
    output:
        allmut = outdir + "/reditools2/{sample}_{type}.redi.{ori}.{ref}.allmut.txt",
        fourbase = outdir + "/reditools2/{sample}_{type}.redi.{ori}.{ref}.fourbase.txt",
    log:
        outdir + "/reditools2/{sample}.count.log"
    threads: 1
    shell:
        """
        # Extract all mismatch counts
        python redi.extract_mismatch_count.gz.py {input} 1> {output.allmut} 2> {log}
        
        # Extract ACGT base counts
        python redi.extract_base_count.gz.py {input} 1> {output.fourbase} 2>> {log}
        """

################################################################################
# Step 9: Extract SNPs from input control samples
################################################################################

rule my_step_09_snp:
    input:
        outdir + "/reditools2/{insample}_{type}.redi.{ori}.{ref}.txt.gz",
    output:
        outdir + "/snp/{insample}_{type}.redi.{ori}.{ref}.snp.txt",
    log:
        outdir + "/snp/{insample}.snp.log"
    threads: 1
    shell:
        """
        python redi.extract_newsnp.gz.py {input} {output} 2> {log}
        """

################################################################################
# Step 10: Extract all mutations from IP samples
################################################################################

rule my_step_10_mut:
    input:
        outdir + "/reditools2/{ipsample}_{type}.redi.{ori}.{ref}.txt.gz",
    output:
        outdir + "/mut/{ipsample}_{type}.redi.{ori}.{ref}.allmut.txt",
    log:
        outdir + "/mut/{ipsample}.allmut.log"
    threads: 1
    shell:
        """
        python redi.extract_allmut.gz.py {input} {output} 2> {log}
        """

################################################################################
# Step 11: Merge SNP sites across all chromosomes for each input sample
################################################################################

rule my_step_11_extract_1:
    input:
        expand(outdir + "/snp/{{insample}}_{{type}}.redi.{{ori}}.{ref}.snp.txt", ref=refs)
    output:
        site = outdir + "/snp/{insample}_{type}.redi.{ori}.newsnp.site.txt"
    params:
        txt = outdir + "/snp/{insample}_{type}.redi.{ori}.allchr.snp.txt",
    threads: 1
    shell:
        """
        # Concatenate all chromosome SNPs
        cat {input} > {params.txt}
        
        # Extract unique chromosome-position pairs
        cut -f 1,2 {params.txt} > {output.site}
        """

################################################################################
# Step 12: Create unified SNP site list across all input samples
################################################################################

rule my_step_12_extract_2:
    input:
        expand(outdir + "/snp/{insample}_{type}.redi.{{ori}}.newsnp.site.txt", insample=insamples, type=type)
    output:
        outdir + "/snp/redi.{ori}.allsnp.site.txt",
    log:
        outdir + "/snp/redi.{ori}.allsnp.log",
    threads: 1
    shell:
        """
        # Merge and deduplicate all SNP sites
        cat {input} | uniq | sort | uniq 1> {output} 2> {log}
        """

################################################################################
# Step 13: Remove SNPs from IP sample mutations
################################################################################

rule my_step_13_extract_3:
    input:
        all = expand(outdir + "/mut/{{ipsample}}_{{type}}.redi.{{ori}}.{ref}.allmut.txt", ref=refs),
        snp = outdir + "/snp/redi.{ori}.allsnp.site.txt",
    output:
        outdir + "/mut/{ipsample}_{type}.redi.{ori}.mut.txt",
    log:
        outdir + "/mut/{ipsample}_{type}.redi.{ori}.mut.log",
    threads: 1
    shell:
        """
        # Remove known SNPs from mutation list
        cat {input.all} | perl rmnewsnp.pl - {input.snp} 1> {output} 2> {log}
        """

################################################################################
# Step 14: Annotate mutations with methylation context (CpG/CHG/CHH)
################################################################################

rule my_step_14_extract_4:
    input:
        expand(outdir + "/mut/{{ipsample}}_{{type}}.redi.{ori}.mut.txt", ori=ori),
    output:
        outdir + "/mod/{ipsample}_{type}.redi.mod.anno.txt",
    params:
        cpg = refdir + "/ucsc.GRCh38.p13.analysis_set.CpG.site.txt",
        chg = refdir + "/ucsc.GRCh38.p13.analysis_set.CHG.site.txt",
        mut = outdir + "/mut/{ipsample}_{type}.redi.mut.anno.txt",
        count = outdir + "/mod/{ipsample}_{type}.redi.mod.anno.count.txt",
    log:
        outdir + "/mut/{ipsample}_{type}.redi.mut.anno.log",
    threads: 1
    shell:
        """
        # Annotate cytosine context (CpG/CHG/CHH)
        perl split_CpGCHGCHH.pl {params.cpg} {params.chg} {input} 1> {params.mut} 2> {log}
        
        # Extract modification sites
        perl extract_mod_from_mut.pl {params.mut} {wildcards.type} 1> {output} 2>> {log}
        
        # Count modifications by context
        perl count_CpGCHGCHH.pl {output} {wildcards.type} 1> {params.count} 2>> {log}
        """

################################################################################
# End of pipeline
################################################################################
