#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * NGS Pipeline Atlas - reproducible bacterial WGS variant-calling pipeline.
 * Stages: FETCH -> SANITIZE (read-pair integrity filter; raw archival
 *         files carry byte corruption that strict parsers reject) ->
 *         FASTQC -> FASTP (trim/QC) -> ALIGN (BWA) ->
 *         GATE (per-sample QC gates) -> CALL (bcftools, passed samples only) ->
 *         REPORT (custom Python summary; no MultiQC)
 */

params.samples = 'samples.tsv'
params.ref     = 'ref/NC_000913.3.fasta'
params.outdir  = 'results/run1'

process PREPARE_REF {
    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta}*"), emit: ref_index

    script:
    """
    bwa index ${fasta}
    samtools faidx ${fasta}
    """
}

process FETCH {
    publishDir "${params.outdir}/fastq", mode: 'copy'

    tag "$sample"
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample), val(r1_url), val(r2_url), val(r1_md5), val(r2_md5)

    output:
    tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz"), emit: reads

    script:
    """
    # The archive rate-limits sustained fetching (HTTP 403 pages) and
    # transfers occasionally truncate, so download until the MD5 of the
    # payload matches the archive's published MD5, with backoff.
    fetch_checked () {
        url=\$1; md5=\$2; out=\$3
        delay=20
        for attempt in 1 2 3 4 5 6; do
            # '|| true': a hard network failure (reset/timeout) must not
            # abort the script under bash -e; the MD5 check below is the
            # only success criterion. --max-time bounds a trickling
            # connection so the loop, not the proxy, controls retries.
            curl -sL --max-time 280 -o \$out \$url || true
            if echo "\$md5  \$out" | md5sum -c - >/dev/null 2>&1; then
                echo "MD5 verified: \$out"
                return 0
            fi
            echo "fetch attempt \$attempt for \$out failed MD5 check; backing off \${delay}s"
            sleep \$delay
            delay=\$(( delay * 2 ))
            [ \$delay -gt 240 ] && delay=240
        done
        echo "ERROR: could not fetch an MD5-valid copy of \$url" >&2
        return 1
    }
    fetch_checked "${r1_url}" "${r1_md5}" ${sample}_1.fastq.gz
    sleep 5
    fetch_checked "${r2_url}" "${r2_md5}" ${sample}_2.fastq.gz
    """
}

process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    tag "$sample"

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    path("*_fastqc.*"), emit: reports

    script:
    """
    fastqc -t ${task.cpus} --outdir . ${r1} ${r2}
    """
}

process SANITIZE {
    publishDir "${params.outdir}/sanitize", mode: 'copy', pattern: "*.sanitize.tsv"

    tag "$sample"

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}.clean_1.fastq.gz"), path("${sample}.clean_2.fastq.gz"), emit: clean
    path("${sample}.sanitize.tsv"), emit: stats

    script:
    """
    sanitize_pairs.py --sample ${sample} \\
                      --r1 ${r1} --r2 ${r2} \\
                      --out1 ${sample}.clean_1.fastq.gz --out2 ${sample}.clean_2.fastq.gz \\
                      --stats ${sample}.sanitize.tsv
    """
}

process FASTP {
    publishDir "${params.outdir}/fastp", mode: 'copy'

    tag "$sample"

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}.trim_1.fastq.gz"), path("${sample}.trim_2.fastq.gz"), emit: trimmed
    tuple val(sample), path("${sample}.fastp.json"), emit: json
    path("${sample}.fastp.html"), emit: html

    script:
    """
    fastp -i ${r1} -I ${r2} \\
          -o ${sample}.trim_1.fastq.gz -O ${sample}.trim_2.fastq.gz \\
          -w ${task.cpus} --detect_adapter_for_pe \\
          -j ${sample}.fastp.json -h ${sample}.fastp.html
    """
}

process ALIGN {
    publishDir "${params.outdir}/bam", mode: 'copy'

    tag "$sample"

    input:
    tuple val(sample), path(r1), path(r2)
    tuple path(ref), path(ref_idx)

    output:
    tuple val(sample), path("${sample}.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} -R '@RG\\tID:${sample}\\tSM:${sample}' ${ref} ${r1} ${r2} \\
        | samtools sort -@ 1 -o ${sample}.bam -
    samtools index ${sample}.bam
    """
}

process GATE {
    publishDir "${params.outdir}/gate", mode: 'copy', pattern: "*.gate.*"

    tag "$sample"

    input:
    tuple val(sample), path(bam), path(fastp_json)

    output:
    tuple val(sample), path("${sample}.gate.tsv"), path("${sample}.gate.log"), emit: gates
    tuple val(sample), path("${sample}.bam"), emit: passed, optional: true

    script:
    """
    samtools index ${bam}
    samtools flagstat ${bam} > ${sample}.flagstat.txt
    samtools coverage ${bam} > ${sample}.coverage.txt
    gate_qc.py --sample ${sample} \\
               --fastp ${fastp_json} \\
               --flagstat ${sample}.flagstat.txt \\
               --coverage ${sample}.coverage.txt \\
               --out ${sample}.gate.tsv \\
               --log ${sample}.gate.log
    # A STOPPED sample's BAM is removed so it cannot flow into variant calling.
    if grep -qP "^${sample}\\tSTOPPED" ${sample}.gate.tsv; then
        rm ${sample}.bam
    fi
    """
}

process CALL {
    publishDir "${params.outdir}/vcf", mode: 'copy'

    tag "$sample"

    input:
    tuple val(sample), path(bam)
    tuple path(ref), path(ref_idx)

    output:
    tuple val(sample), path("${sample}.vcf.gz"), emit: vcf
    path("${sample}.calls.tsv"), emit: calls

    script:
    """
    samtools index ${bam}
    bcftools mpileup -f ${ref} -Ou -q 20 -Q 20 --max-depth 5000 ${bam} \\
        | bcftools call -mv -Oz -o ${sample}.vcf.gz
    bcftools index -t ${sample}.vcf.gz
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\n' ${sample}.vcf.gz \\
        | sort -k1,1 -k2,2n > ${sample}.calls.tsv
    """
}

process REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path gate_files
    path call_files

    output:
    path("cohort_report.txt"), emit: report
    path("summary.tsv"), emit: summary
    path("gate_log_all.txt"), emit: gate_log

    script:
    """
    mkdir -p gates calls
    for f in ${gate_files}; do case "\$f" in *.gate.tsv) cp \$f gates/ ;; *.gate.log) cp \$f gates/ ;; esac; done
    for f in ${call_files}; do cp \$f calls/; done
    cat gates/*.gate.log > gate_log_all.txt
    make_report.py --gate-dir gates --calls-dir calls \\
                   --out-report cohort_report.txt --out-summary summary.tsv
    """
}

workflow {
    ref_ch = PREPARE_REF(file(params.ref))

    samples_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample, row.r1_url, row.r2_url, row.r1_md5, row.r2_md5) }

    fetched = FETCH(samples_ch)
    sanitized = SANITIZE(fetched.reads)
    FASTQC(sanitized.clean)
    fastp_out = FASTP(sanitized.clean)
    aligned = ALIGN(fastp_out.trimmed, ref_ch.ref_index)
    gate_in = aligned.bam.join(fastp_out.json)
    gate_out = GATE(gate_in)
    called = CALL(gate_out.passed, ref_ch.ref_index)
    REPORT(gate_out.gates.map { s, tsv, log -> [tsv, log] }.flatten().collect(),
           called.calls.collect())
}
