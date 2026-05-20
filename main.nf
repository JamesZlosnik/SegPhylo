#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================
// Parameters
// ============================================================

params.l_fasta        = null
params.m_fasta        = null
params.s_fasta        = null
params.root_l         = null
params.root_m         = null
params.root_s         = null
params.outdir         = "results"
params.mafft_args      = "--auto"
params.iqtree_model    = "GTR+G"
params.iqtree_boot     = 1000
params.label_taxa_max  = 250    // suppress sample labels in barcode plot above this count
params.min_len_l       = 6000   // minimum length for segment L
params.min_len_m       = 3000   // minimum length for segment M
params.min_len_s       = 1000   // minimum length for segment S
params.mask_from_start = 30     // positions to mask from start of each aligned segment
params.mask_from_end   = 50     // positions to mask from end of each aligned segment
params.help            = false

// ============================================================
// Help
// ============================================================

def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║           Segmented Virus Phylogenetics Pipeline         ║
    ╚══════════════════════════════════════════════════════════╝

    Builds per-segment and concatenated phylogenetic trees from
    a tri-segmented virus (L, M, S). Headers must follow the
    format:  >SampleName|SEGMENT

    Usage:
        nextflow run main.nf [options]

    Required:
        --l_fasta       Multifasta for segment L
        --m_fasta       Multifasta for segment M
        --s_fasta       Multifasta for segment S
        --root_l        Root/outgroup sequence for segment L
        --root_m        Root/outgroup sequence for segment M
        --root_s        Root/outgroup sequence for segment S

    Optional:
        --outdir            Output directory            [default: results]
        --mafft_args        MAFFT alignment arguments  [default: --auto]
        --iqtree_model      Substitution model for per-segment trees
                                                        [default: GTR+G]
        --iqtree_boot       Ultrafast bootstrap reps   [default: 1000]
        --min_len_l         Minimum length for segment L [default: 6000]
        --min_len_m         Minimum length for segment M [default: 3000]
        --min_len_s         Minimum length for segment S [default: 1000]
        --mask_from_start   Positions to mask from start of each aligned segment [default: 30]
        --mask_from_end     Positions to mask from end of each aligned segment   [default: 50]
        --label_taxa_max    Hide sample labels in barcode plot above
                            this number of taxa         [default: 250]

    Profiles:
        -profile conda          Auto-build conda environment (no defaults channel)
        -profile slurm          Submit jobs via SLURM
        -profile slurm,conda    Combine profiles

    Outputs:
        results/01_filtered/        Filtered fastas + filtering report
        results/02_alignments/      Per-segment and concatenated alignments
        results/03_trees/           Concatenated tree (primary output) + per-segment trees
        results/04_visualization/   Barcode alignment plot (PNG + SVG)

    Example:
        nextflow run main.nf \\
            --l_fasta seqs_L.fasta --m_fasta seqs_M.fasta --s_fasta seqs_S.fasta \\
            --root_l root_L.fasta  --root_m root_M.fasta  --root_s root_S.fasta  \\
            --outdir results -profile conda
    """.stripIndent()
}

// ============================================================
// Processes
// ============================================================

/*
 * Filter: retain only samples present in all three segments (L, M, S).
 * Expects headers in the format:  >SampleName|SEGMENT
 * Produces a filtering report listing retained and discarded samples.
 */
process FILTER_COMPLETE_SAMPLES {
    tag "filter"
    publishDir "${params.outdir}/01_filtered", mode: 'copy'

    input:
    path l_fasta
    path m_fasta
    path s_fasta

    output:
    path "filtered_L.fasta",     emit: l
    path "filtered_M.fasta",     emit: m
    path "filtered_S.fasta",     emit: s
    path "filtering_report.txt", emit: report

    script:
    """
    filter_complete_samples.py \\
        --l_fasta   ${l_fasta} \\
        --m_fasta   ${m_fasta} \\
        --s_fasta   ${s_fasta} \\
        --min_len_l ${params.min_len_l} \\
        --min_len_m ${params.min_len_m} \\
        --min_len_s ${params.min_len_s}
    """
}

/*
 * Align each segment against the root sequence using reference-guided MAFFT
 * (--addfragments). The root acts as the fixed reference; samples are aligned
 * to it rather than being pooled for de novo alignment. --keeplength prevents
 * extra columns being inserted relative to the reference.
 * The root (OUTGROUP) appears first in the output alignment.
 */
process MAFFT_ALIGN {
    tag "align: ${segment}"
    publishDir "${params.outdir}/02_alignments", mode: 'copy', pattern: "*_aligned.fasta"

    input:
    tuple val(segment), path(fasta), path(root_fasta)

    output:
    tuple val(segment), path("${segment}_aligned.fasta")

    script:
    """
    # Write outgroup-only reference fasta (handles FASTA or GenBank input)
    add_outgroup.py \\
        --segment      ${segment} \\
        --fasta        ${fasta} \\
        --root_fasta   ${root_fasta} \\
        --outgroup_only

    # Reference-guided alignment: add samples to the single outgroup reference.
    # --keeplength keeps alignment width equal to the reference length.
    mafft ${params.mafft_args} --keeplength --thread ${task.cpus} \\
        --addfragments ${fasta} outgroup_ref.fasta > ${segment}_aligned.fasta
    """
}

/*
 * Mask terminal positions of each aligned segment by replacing with N.
 * Mirrors augur mask defaults: 30 bp from start, 50 bp from end.
 * These regions are often poorly covered and can distort tree topology.
 */
process MASK_ALIGNMENT {
    tag "mask: ${segment}"
    publishDir "${params.outdir}/02_alignments", mode: 'copy'

    input:
    tuple val(segment), path(alignment)

    output:
    tuple val(segment), path("${segment}_masked.fasta")

    script:
    """
    mask_alignment.py \\
        --alignment       ${alignment} \\
        --segment         ${segment} \\
        --mask_from_start ${params.mask_from_start} \\
        --mask_from_end   ${params.mask_from_end}
    """
}

/*
 * Build an IQ-TREE2 phylogenetic tree for each individual segment.
 * OUTGROUP is used as the root.
 */
process IQTREE_SEGMENT {
    tag "tree: ${segment}"
    publishDir "${params.outdir}/03_trees/per_segment", mode: 'copy'

    input:
    tuple val(segment), path(alignment)

    output:
    tuple val(segment), path("${segment}.nwk"), emit: tree
    path "${segment}.*",                         emit: all_files

    script:
    """
    iqtree2 \\
        -s  ${alignment} \\
        --prefix ${segment} \\
        -o  OUTGROUP \\
        -m  ${params.iqtree_model} \\
        -B  ${params.iqtree_boot} \\
        -czb -st DNA \\
        --redo \\
        -T  ${task.cpus}
    mv ${segment}.treefile ${segment}.nwk
    """
}

/*
 * Concatenate the three per-segment alignments in L->M->S order.
 * Also generates an IQ-TREE2-compatible RAxML-style partition file so that
 * each segment gets its own substitution model in the concatenated tree.
 */
process CONCATENATE_ALIGNMENTS {
    tag "concatenate"
    publishDir "${params.outdir}/02_alignments", mode: 'copy'

    input:
    path l_aln
    path m_aln
    path s_aln

    output:
    path "concatenated_aligned.fasta", emit: fasta
    path "partitions.txt",             emit: partition

    script:
    """
    concatenate_alignments.py \\
        --l_aln ${l_aln} \\
        --m_aln ${m_aln} \\
        --s_aln ${s_aln}
    """
}

/*
 * Build the full concatenated IQ-TREE2 tree using a partition model
 * (one GTR+G substitution model per segment).
 */
process IQTREE_CONCATENATED {
    tag "concatenated tree"
    publishDir "${params.outdir}/03_trees", mode: 'copy'

    input:
    path alignment
    path partition

    output:
    path "concatenated.nwk",  emit: tree
    path "concatenated.*",    emit: all_files

    script:
    """
    iqtree2 \\
        -s  ${alignment} \\
        -p  ${partition} \\
        --prefix concatenated \\
        -o  OUTGROUP \\
        -B  ${params.iqtree_boot} \\
        -czb -st DNA \\
        --redo \\
        -T  ${task.cpus}
    mv concatenated.treefile concatenated.nwk
    """
}

/*
 * Barcode-style visualization of the concatenated alignment.
 * Differences are shown relative to the OUTGROUP sequence.
 * Outputs: alignment_barcode.png + alignment_barcode.svg
 */
process VISUALIZE_ALIGNMENT {
    tag "visualize"
    publishDir "${params.outdir}/04_visualization", mode: 'copy'

    input:
    path alignment
    path partition

    output:
    path "alignment_barcode.png", emit: png
    path "alignment_barcode.svg", emit: svg

    script:
    """
    visualize_alignment.py \\
        --alignment       ${alignment} \\
        --partition       ${partition} \\
        --mask_from_start ${params.mask_from_start} \\
        --mask_from_end   ${params.mask_from_end} \\
        --label_taxa_max  ${params.label_taxa_max}
    """
}

// ============================================================
// Workflow
// ============================================================

workflow {

    // --- Help / validation ---
    def required = ['l_fasta', 'm_fasta', 's_fasta', 'root_l', 'root_m', 'root_s']

    // Show help if requested, or if run with no arguments at all
    if (params.help || required.every { !params[it] }) {
        helpMessage()
        exit 0
    }

    // Report all missing required params at once rather than stopping at the first
    def missing = required.findAll { !params[it] }
    if (missing) {
        log.error "Missing required parameter(s): ${missing.collect { "--${it}" }.join(', ')}"
        log.info  "Run with --help for usage information."
        exit 1
    }

    // --- Inputs ---
    l_fasta = Channel.fromPath(params.l_fasta, checkIfExists: true)
    m_fasta = Channel.fromPath(params.m_fasta, checkIfExists: true)
    s_fasta = Channel.fromPath(params.s_fasta, checkIfExists: true)
    root_l  = Channel.fromPath(params.root_l,  checkIfExists: true)
    root_m  = Channel.fromPath(params.root_m,  checkIfExists: true)
    root_s  = Channel.fromPath(params.root_s,  checkIfExists: true)

    // --- Step 1: Filter to complete trios ---
    FILTER_COMPLETE_SAMPLES(l_fasta, m_fasta, s_fasta)

    // --- Step 2: Pair each filtered segment with its root sequence ---
    // Channels emit: [ segment_label, filtered_fasta ]
    seg_filtered_ch = FILTER_COMPLETE_SAMPLES.out.l.map { ["L", it] }
        .mix( FILTER_COMPLETE_SAMPLES.out.m.map { ["M", it] } )
        .mix( FILTER_COMPLETE_SAMPLES.out.s.map { ["S", it] } )

    // Channels emit: [ segment_label, root_fasta ]
    root_ch = root_l.map { ["L", it] }
        .mix( root_m.map { ["M", it] } )
        .mix( root_s.map { ["S", it] } )

    // Join on segment label -> [ segment_label, filtered_fasta, root_fasta ]
    seg_root_ch = seg_filtered_ch.join(root_ch)

    // --- Step 3: Reference-guided alignment then mask terminals ---
    MAFFT_ALIGN(seg_root_ch)
    MASK_ALIGNMENT(MAFFT_ALIGN.out)

    // --- Step 4: Per-segment trees (on masked alignments) ---
    IQTREE_SEGMENT(MASK_ALIGNMENT.out)

    // --- Step 5: Concatenate masked alignments ---
    l_aln = MASK_ALIGNMENT.out.filter { it[0] == "L" }.map { it[1] }
    m_aln = MASK_ALIGNMENT.out.filter { it[0] == "M" }.map { it[1] }
    s_aln = MASK_ALIGNMENT.out.filter { it[0] == "S" }.map { it[1] }

    CONCATENATE_ALIGNMENTS(l_aln, m_aln, s_aln)

    // --- Step 6: Concatenated tree ---
    IQTREE_CONCATENATED(
        CONCATENATE_ALIGNMENTS.out.fasta,
        CONCATENATE_ALIGNMENTS.out.partition
    )

    // --- Step 7: Barcode alignment visualization ---
    VISUALIZE_ALIGNMENT(
        CONCATENATE_ALIGNMENTS.out.fasta,
        CONCATENATE_ALIGNMENTS.out.partition
    )
}
