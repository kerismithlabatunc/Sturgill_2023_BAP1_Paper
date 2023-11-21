#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//
// Runs a modified version of the UNC Lineberger somatic DNA workflow
// NOTE: Will not run without valid paths to software and inputs used
//

// Basic input and output information
params.cohort_name = 'cohort'
params.input = ''
params.output_dir = ''
params.systemlog = ''
params.mail = ''

params.publish_dir_mode = 'copy'
publish_dir_mode = params.publish_dir_mode

params.publish_bams = true
publish_bams = params.publish_bams

cohort_name = params.cohort_name
params.species = 'homo_sapiens'
params.vep_cache_version = 103
params.vep_config = ''

//TODO: Need to provision these
params.vep_retain_ann = ''

if (params.output_dir != '') {
    output_dir = params.output_dir
} else {
    println "--output_dir PATH must be specified!"
    System.exit(-1)
}

// set up reference genome

params.ref_genome_base = 'GRCh38'
params.ref_path = '/path/to/GRCh38.d1.vd1'
params.ref_fasta = 'GRCh38.d1.vd1.fa'

ref_genome = params.ref_genome_base
ref_genome_dir = file(params.ref_path)
ref_genome_fasta = params.ref_fasta
ref_genome_dict = ref_genome_fasta.take(ref_genome_fasta.lastIndexOf('.fa')) + ".dict"
ref_genome_fai = ref_genome_fasta + '.fai'

params.ref_bedtools_genome = ref_genome_fasta + '.bedtoolsGenome'
ref_bedtools_genome = params.ref_bedtools_genome

// Program path defaults:

params.abra_jar_path = '/opt/conda/share/abra2-2.24-0'
params.vep_exec_path = '/opt/conda/bin'
vep_exec_path = params.vep_exec_path
params.vep_cache = ''
vep_cache = file (params.vep_cache)

params.bin_path = '/opt/bin'
bin_path = params.bin_path
params.data_path = '/opt/data'
data_path = params.data_path

// default coverage regions:  gencode V32
params.bait_name = 'gencode_v32'
bait_name = params.bait_name
params.capture_bed = '/path/to/gencode.v32.exons.merged.sort.slop10.bed'
capture_bed = file(params.capture_bed)

params.target_bed = '/path/to/gencode.v32.exons.merged.sort.bed'
target_bed = file(params.target_bed)

params.target_bed_bgz = '/path/to/bap1_project_targets.bed.gz'
target_bed_bgz = file(params.target_bed_bgz)

//  quality filter params:
params.cadabra_qual = 13.3
params.indel_strelka2_somatic_evs = 17.3
params.snv_strelka2_somatic_evs = 20.0

// Handle variable input types:

// Initialize input channels
bam_channel_input = Channel.empty()

// set up bam_channel
println "BAM input!"
bam_channel_input = Channel
        .fromPath(file(params.input))
        .splitCsv(header:true, sep:'\t')
        .map{ row-> tuple(row.subject, file(row.tumor_bam), file(row.normal_bam))}

/////////////////////////////////////////////////////////
// Per cohort prep processes.
/////////////////////////////////////////////////////////

// Index

process samtools_index {
    cpus 1
    memory '2 GB'
    errorStrategy {task.attempt < 4 ? 'retry': 'finish'}
    maxRetries 4

    tag "${cohort_name}+${subject}"

    cache 'lenient'

    input:
        tuple val(subject),
            path(tumor_bam),
            path(normal_bam)

    output:
        tuple val(subject),
            path("${tumor_bam}.bai"), emit: tumor_bai
        tuple val(subject),
            path("${normal_bam}.bai"), emit: normal_bai

    """
    samtools index ${tumor_bam}
    samtools index ${normal_bam}
    """

}

// use ABRA2 to realign the Tumor and Normal for downstream somatic variant calling.

process abra {
    cpus 8
    memory {18.GB * task.attempt}
    errorStrategy {task.attempt < 4 ? 'retry': 'finish'}
    maxRetries 4

    tag "${cohort_name}+${subject}"

    cache 'lenient'

    input:
        path ref_genome_dir
        val ref_genome_fasta
        path target_bed
        tuple val(subject),
            path(tumor_bam),
            path(normal_bam)

    output:
        tuple val(subject),
            path('tumor.abra.bam'), emit: tumor_abra_bam_sortindex
        tuple val(subject),
            path('normal.abra.bam'), emit: normal_abra_bam_sortindex

    script:
        javaMem = task.attempt * 15

    """
    java -Xmx${javaMem}g -Xms${javaMem}g \
        -jar ${params.abra_jar_path}/abra2.jar \
        --threads 8 \
        --tmpdir . \
        --out tumor.abra.bam,normal.abra.bam \
        --in ${tumor_bam},${normal_bam} \
        --ref ${ref_genome_dir}/${ref_genome_fasta} \
        --undup \
        --no-edge-ci \
        --targets ${target_bed}
    """
}

// Mark duplicates on the abrafied BAM files for T/N separately.

process mark_dupes_abra_tumor {
    cpus 16
    memory '32 GB'
    tag "${cohort_name}+${subject}_T"

    publishDir "${output_dir}/${subject}",
        mode: publish_dir_mode,
        pattern: '*.{bam,bai,txt}',
        enabled: publish_bams

    input:
        tuple val(subject),
            path(tumor_bam)
    output:
        tuple val(subject),
            path('tumor.abra.md.bam'),
            path('tumor.abra.md.bam.bai'), emit: tumor_abra_md_bam
        path '*.metrics.txt', emit: tumor_abra_metrics

    script:

    bamZipLevel = publish_bams ? "9" : "1"

    """
    bammarkduplicates \
        I=${tumor_bam} \
        O=tumor.abra.md.bam \
        markthreads=${task.cpus} \
        M=tumor.abra.md.metrics.txt \
        index=1 \
        level=${bamZipLevel} \
        rmdup=0
    """
}

process mark_dupes_abra_normal {
    cpus 16
    memory '32 GB'
    tag "${cohort_name}+${subject}_N"

    publishDir "${output_dir}/${subject}",
        mode: publish_dir_mode,
        pattern: '*.{bam,bai,txt}',
        enabled: publish_bams

    input:
        tuple val(subject),
            path(normal_bam)
    output:
        tuple val(subject),
            path('normal.abra.md.bam'),
            path('normal.abra.md.bam.bai'), emit: normal_abra_md_bam
        path '*.metrics.txt', emit: normal_abra_metrics

    script:

    bamZipLevel = publish_bams ? "9" : "1"

    """
    bammarkduplicates \
        I=${normal_bam} \
        O=normal.abra.md.bam \
        markthreads=${task.cpus} \
        M=normal.abra.md.metrics.txt \
        index=1 \
        level=${bamZipLevel} \
        rmdup=0
    """
}


// call variants

process cadabra {
    cpus 8
    memory {18.GB * task.attempt}
    errorStrategy {task.attempt < 4 ? 'retry': 'finish'}
    maxRetries 4

    tag "${cohort_name}+${subject}"

    input:
        path ref_genome_dir
        val ref_genome_fasta
        tuple val(subject),
            path(tumor_bam),
            path(tumor_bai),
            path(normal_bam),
            path(normal_bai)

    output:
        tuple val(subject),
            path('cadabra.vcf'), emit: cadabra_calls

    script:
        javaMem = task.attempt * 16

    """
    java -Xmx${javaMem}g -XX:+UseParallelGC \
        -cp ${params.abra_jar_path}/abra2.jar abra.cadabra.Cadabra \
        --threads ${task.cpus} \
        --tumor ${tumor_bam} \
        --normal ${normal_bam} \
        --ref ${ref_genome_dir}/${ref_genome_fasta} > cadabra.vcf
    """
}

process strelka_somatic {
    cpus 8
    memory '8 GB'
    tag "${cohort_name}+${subject}"

    input:
        path ref_genome_dir
        val ref_genome_fasta
        tuple val(subject),
            path(tumor_bam),
            path(tumor_bai),
            path(normal_bam),
            path(normal_bai)

    output:
        tuple val(subject),
            path('strelka2.snvs.vcf'), emit: strelka2_snvs
        tuple val(subject),
            path('strelka2.indels.vcf'), emit: strelka2_indels

    """
    configureStrelkaSomaticWorkflow.py \
        --targeted \
        --tumorBam ${tumor_bam} \
        --normalBam ${normal_bam} \
        --runDir StrelkaSomaticWorkflow \
        --referenceFasta ${ref_genome_dir}/${ref_genome_fasta} \
        --callRegions ${target_bed_bgz}

    cd StrelkaSomaticWorkflow

    ./runWorkflow.py \
        -m local \
        -j ${task.cpus}

    cd ..

    rm -rfv StrelkaSomaticWorkflow/workspace/pyflow.data/logs/tmp

    mv -v StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz       strelka2.snvs.vcf.gz
    mv -v StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz     strelka2.indels.vcf.gz
    gzip -d strelka2.snvs.vcf.gz
    gzip -d strelka2.indels.vcf.gz
    """
}

process vcf2maf_cadabra {
    cpus 16
    memory '16 GB'
    tag "${cohort_name}+${subject}"

    publishDir "${output_dir}/${subject}",
        mode: publish_dir_mode,
        pattern: '*.maf'

    input:
        path ref_genome_dir
        val ref_genome_fasta
        tuple val(subject),
            path(cadabra_vcf)

        path vep_cache
        val vep_config
        val vep_retain_ann

    output:
        tuple val(subject),
            path('cadabra.maf'), emit: cadabra_maf

    script:

    if (ref_genome == 'hg19'){
        ref_genome = 'GRCh37'
    }

    vep_config_cmd = ''

    if (vep_config) {
        vep_config_cmd += "--vep-config " + vep_config
        vep_config_cmd += " --retain-ann " + vep_retain_ann
    }

    """
    set -ueo pipefail

    vcf2maf.pl --vep-path ${vep_exec_path} \
        --vep-data ${vep_cache} \
        --vep-forks ${task.cpus} \
        --ref-fasta ${ref_genome_dir}/${ref_genome_fasta} \
        --ncbi-build ${ref_genome} \
        --input-vcf ${cadabra_vcf} \
        --output-maf cadabra.maf \
        ${vep_config_cmd} \
        --tumor-id TUMOR \
        --normal-id NORMAL \
        --species ${params.species} \
        --cache-version ${params.vep_cache_version}
    """
}

process vcf2maf_strelka2 {
    cpus 16
    memory '16 GB'
    tag "${cohort_name}+${subject}"

    publishDir "${output_dir}/${subject}",
        mode: publish_dir_mode,
        pattern: '*.maf'

    input:
        path ref_genome_dir
        val ref_genome_fasta
        tuple val(subject),
            path(strelka2_snvs),
            path(strelka2_indels)

        path vep_cache
        val vep_config
        val vep_retain_ann

    output:
        tuple val(subject),
            path('strelka2.snvs.maf'),
            path('strelka2.indels.maf'), emit: strelka2_maf

    script:

    if (ref_genome == 'hg19'){
        ref_genome = 'GRCh37'
    }

    vep_config_cmd = ''

    if (vep_config) {
        vep_config_cmd += "--vep-config " + vep_config
        vep_config_cmd += " --retain-ann " + vep_retain_ann
    }

    """
    set -ueo pipefail

    vcf2maf.pl --vep-path ${vep_exec_path} \
        --vep-data ${vep_cache} \
        --vep-forks ${task.cpus} \
        --ref-fasta ${ref_genome_dir}/${ref_genome_fasta} \
        --ncbi-build ${ref_genome} \
        --input-vcf ${strelka2_snvs} \
        --output-maf strelka2.snvs.maf \
        ${vep_config_cmd} \
        --tumor-id TUMOR \
        --normal-id NORMAL \
        --species ${params.species} \
        --cache-version ${params.vep_cache_version}

    vcf2maf.pl --vep-path ${vep_exec_path} \
        --vep-data ${vep_cache} \
        --vep-forks ${task.cpus} \
        --ref-fasta ${ref_genome_dir}/${ref_genome_fasta} \
        --ncbi-build ${ref_genome} \
        --input-vcf ${strelka2_indels} \
        --output-maf strelka2.indels.maf \
        ${vep_config_cmd} \
        --tumor-id TUMOR \
        --normal-id NORMAL \
        --species ${params.species} \
        --cache-version ${params.vep_cache_version}
    """
}

workflow.onComplete {
    def msg = """\
        ${workflow.manifest.name} workflow execution summary
        ---------------------------
        Project : ${workflow.projectDir}
        Script  : ${workflow.scriptFile}
        Git info: ${workflow.repository} - ${workflow.revision} [${workflow.commitId}]
        Cmd line: ${workflow.commandLine}
        Manifest's pipeline version: ${workflow.manifest.version}
           --
        Started at  : ${workflow.start}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        Resumed?    : ${workflow.resume}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        ---------------------------
        """
        .stripIndent()

    if (params.mail !=null && params.mail.length()>0) {
        println ("Sending report to ${params.mail}")
        sendMail(to: "${params.mail}",
            subject: "My ${workflow.manifest.name} execution",
            body: msg)
    }

    if (params.systemlog != null && params.systemlog.length()>0) {
        println ("Logging run info to ${params.systemlog}")
        systemlog = file (params.systemlog)
        systemlog.append(msg)
    }
    println (msg)
}

println "Start workflow"
workflow {
    samtools_index(bam_channel_input)
    abra(ref_genome_dir,ref_genome_fasta,target_bed,bam_channel_input)
    mark_dupes_abra_tumor(abra.out.tumor_abra_bam_sortindex)
    mark_dupes_abra_normal(abra.out.normal_abra_bam_sortindex)
    strelka_somatic(ref_genome_dir,ref_genome_fasta,mark_dupes_abra_tumor.out.tumor_abra_md_bam.join(mark_dupes_abra_normal.out.normal_abra_md_bam))
    cadabra(ref_genome_dir,ref_genome_fasta,mark_dupes_abra_tumor.out.tumor_abra_md_bam.join(mark_dupes_abra_normal.out.normal_abra_md_bam))
    vcf2maf_cadabra(ref_genome_dir,ref_genome_fasta,cadabra.out.cadabra_calls,vep_cache,params.vep_config,params.vep_retain_ann)
    vcf2maf_strelka2(ref_genome_dir,ref_genome_fasta,strelka_somatic.out.strelka2_snvs.join(strelka_somatic.out.strelka2_indels),vep_cache,params.vep_config,params.vep_retain_ann)
}


