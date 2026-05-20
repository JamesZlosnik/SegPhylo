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
params.dropped_strains = null   // optional file of strain names to exclude (one per line)
params.mask_from_start = 30     // positions to mask from start of each aligned segment
params.mask_from_end   = 50     // positions to mask from end of each aligned segment
params.remove_reference = false // remove OUTGROUP from alignment before tree building
                                // (mirrors augur align --remove-reference)
params.root            = "outgroup" // tree rooting: "outgroup" or "midpoint"
params.segment_trees_all_passing = false // build per-segment trees from all
                                         // length-passing samples for that segment
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
        --iqtree_model      Substitution model for per-segment and concatenated trees
                                                        [default: GTR+G]
        --iqtree_boot       Ultrafast bootstrap reps   [default: 1000]
        --min_len_l         Minimum length for segment L [default: 6000]
        --min_len_m         Minimum length for segment M [default: 3000]
        --min_len_s         Minimum length for segment S [default: 1000]
        --dropped_strains   Optional file of sample names to exclude, one per line
                            (comment lines starting with '#' are ignored)
        --mask_from_start   Positions to mask from start of each aligned segment [default: 30]
        --mask_from_end     Positions to mask from end of each aligned segment   [default: 50]
        --remove_reference  Remove OUTGROUP from alignment before tree building  [default: false]
                            Mirrors augur align --remove-reference
        --root              Tree rooting method: 'outgroup' or 'midpoint'        [default: outgroup]
                            'midpoint' mirrors augur refine --root mid_point
                            Note: --root outgroup requires --remove_reference false
        --segment_trees_all_passing
                            Build per-segment trees from all samples passing that
                            segment's threshold, while concatenated tree still
                            uses complete trios only                  [default: false]
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
 * Also emits per-segment files containing all samples that pass the threshold
 * for that segment, for optional per-segment tree building.
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
    path dropped_strains

    output:
    path "filtered_L.fasta",     emit: l
    path "filtered_M.fasta",     emit: m
    path "filtered_S.fasta",     emit: s
    path "segment_filtered_L.fasta", emit: segment_l
    path "segment_filtered_M.fasta", emit: segment_m
    path "segment_filtered_S.fasta", emit: segment_s
    path "filtering_report.txt", emit: report

    script:
    def drop_arg = dropped_strains.name != 'NO_FILE' ? "--dropped_strains ${dropped_strains}" : ""
    """
    python ${projectDir}/bin/filter_complete_samples.py \\
        --l_fasta   ${l_fasta} \\
        --m_fasta   ${m_fasta} \\
        --s_fasta   ${s_fasta} \\
        --min_len_l ${params.min_len_l} \\
        --min_len_m ${params.min_len_m} \\
        --min_len_s ${params.min_len_s} \\
        ${drop_arg}
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
    python ${projectDir}/bin/add_outgroup.py \\
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
    python ${projectDir}/bin/mask_alignment.py \\
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
    def outgroup_arg = params.root == "outgroup" ? "-o OUTGROUP" : ""
    """
    iqtree2 \\
        -s  ${alignment} \\
        --prefix ${segment} \\
        ${outgroup_arg} \\
        -m  ${params.iqtree_model} \\
        -B  ${params.iqtree_boot} \\
        -czb -st DNA \\
        --redo \\
        -T  ${task.cpus}
    mv ${segment}.treefile ${segment}.nwk
    """
}

/*
 * Remove the OUTGROUP sequence from an alignment before tree building.
 * Mirrors augur align --remove-reference. Only invoked when --remove_reference true.
 */
process STRIP_OUTGROUP {
    tag "strip outgroup: ${segment}"

    input:
    tuple val(segment), path(alignment)

    output:
    tuple val(segment), path("${segment}_stripped.fasta")

    script:
    """
    python ${projectDir}/bin/strip_outgroup.py \\
        --alignment ${alignment} \\
        --output    ${segment}_stripped.fasta
    """
}

/*
 * Remove the OUTGROUP sequence from the concatenated alignment before tree building.
 * Only invoked when --remove_reference true.
 */
process STRIP_OUTGROUP_CONCAT {
    tag "strip outgroup: concatenated"

    input:
    path alignment

    output:
    path "concatenated_stripped.fasta"

    script:
    """
    python ${projectDir}/bin/strip_outgroup.py \\
        --alignment ${alignment} \\
        --output    concatenated_stripped.fasta
    """
}

/*
 * Apply midpoint rooting to a per-segment tree.
 * Only invoked when --root midpoint.
 */
process MIDPOINT_ROOT {
    tag "midpoint root: ${segment}"
    publishDir "${params.outdir}/03_trees/per_segment", mode: 'copy'

    input:
    tuple val(segment), path(tree)

    output:
    tuple val(segment), path("${segment}.nwk")

    script:
    """
    python ${projectDir}/bin/midpoint_root.py \\
        --tree   ${tree} \\
        --output ${segment}.nwk
    """
}

/*
 * Apply midpoint rooting to the concatenated tree.
 * Only invoked when --root midpoint.
 */
process MIDPOINT_ROOT_CONCAT {
    tag "midpoint root: concatenated"
    publishDir "${params.outdir}/03_trees", mode: 'copy'

    input:
    path tree

    output:
    path "concatenated.nwk"

    script:
    """
    python ${projectDir}/bin/midpoint_root.py \\
        --tree   ${tree} \\
        --output concatenated.nwk
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
    python ${projectDir}/bin/concatenate_alignments.py \\
        --l_aln ${l_aln} \\
        --m_aln ${m_aln} \\
        --s_aln ${s_aln} \\
        --model ${params.iqtree_model}
    """
}

/*
 * Build the full concatenated IQ-TREE2 tree using a partition model
 * (one substitution model per segment).
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
    def outgroup_arg = params.root == "outgroup" ? "-o OUTGROUP" : ""
    """
    iqtree2 \\
        -s  ${alignment} \\
        -p  ${partition} \\
        --prefix concatenated \\
        ${outgroup_arg} \\
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
    python ${projectDir}/bin/visualize_alignment.py \\
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

    // Validate root / remove_reference combination
    if (params.root == "outgroup" && params.remove_reference) {
        log.error "--root outgroup cannot be used with --remove_reference true: " +
                  "the OUTGROUP sequence will not be in the alignment for rooting. " +
                  "Use --root midpoint when removing the reference."
        exit 1
    }
    if (!["outgroup", "midpoint"].contains(params.root)) {
        log.error "--root must be 'outgroup' or 'midpoint' (got: '${params.root}')"
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
    dropped_strains = params.dropped_strains
        ? Channel.fromPath(params.dropped_strains, checkIfExists: true)
        : Channel.fromPath("NO_FILE", checkIfExists: false).ifEmpty { file("NO_FILE") }

    FILTER_COMPLETE_SAMPLES(l_fasta, m_fasta, s_fasta, dropped_strains)

    // --- Step 2: Pair each filtered segment with its root sequence ---
    // Channels emit: [ segment_label, filtered_fasta ]
    complete_seg_ch = FILTER_COMPLETE_SAMPLES.out.l.map { ["L", it] }
        .mix( FILTER_COMPLETE_SAMPLES.out.m.map { ["M", it] } )
        .mix( FILTER_COMPLETE_SAMPLES.out.s.map { ["S", it] } )

    all_passing_seg_ch = FILTER_COMPLETE_SAMPLES.out.segment_l.map { ["L", it] }
        .mix( FILTER_COMPLETE_SAMPLES.out.segment_m.map { ["M", it] } )
        .mix( FILTER_COMPLETE_SAMPLES.out.segment_s.map { ["S", it] } )

    seg_filtered_ch = params.segment_trees_all_passing ? all_passing_seg_ch : complete_seg_ch

    // Channels emit: [ segment_label, root_fasta ]
    root_ch = root_l.map { ["L", it] }
        .mix( root_m.map { ["M", it] } )
        .mix( root_s.map { ["S", it] } )

    // Join on segment label -> [ segment_label, filtered_fasta, root_fasta ]
    seg_root_ch = seg_filtered_ch.join(root_ch)

    // --- Step 3: Reference-guided alignment then mask terminals ---
    MAFFT_ALIGN(seg_root_ch)
    MASK_ALIGNMENT(MAFFT_ALIGN.out)

    // --- Step 4: Optionally strip OUTGROUP, then build per-segment trees ---
    def seg_tree_input_ch
    if (params.remove_reference) {
        STRIP_OUTGROUP(MASK_ALIGNMENT.out)
        seg_tree_input_ch = STRIP_OUTGROUP.out
    } else {
        seg_tree_input_ch = MASK_ALIGNMENT.out
    }

    IQTREE_SEGMENT(seg_tree_input_ch)

    if (params.root == "midpoint") {
        MIDPOINT_ROOT(IQTREE_SEGMENT.out.tree)
    }

    // --- Step 5: Concatenate masked alignments ---
    // Always concatenate the full (with-OUTGROUP) alignments for visualization.
    // Tree building uses the stripped version if --remove_reference is set.
    l_aln = MASK_ALIGNMENT.out.filter { it[0] == "L" }.map { it[1] }
    m_aln = MASK_ALIGNMENT.out.filter { it[0] == "M" }.map { it[1] }
    s_aln = MASK_ALIGNMENT.out.filter { it[0] == "S" }.map { it[1] }

    CONCATENATE_ALIGNMENTS(l_aln, m_aln, s_aln)

    // --- Step 6: Concatenated tree ---
    def concat_tree_input
    if (params.remove_reference) {
        STRIP_OUTGROUP_CONCAT(CONCATENATE_ALIGNMENTS.out.fasta)
        concat_tree_input = STRIP_OUTGROUP_CONCAT.out
    } else {
        concat_tree_input = CONCATENATE_ALIGNMENTS.out.fasta
    }

    IQTREE_CONCATENATED(concat_tree_input, CONCATENATE_ALIGNMENTS.out.partition)

    if (params.root == "midpoint") {
        MIDPOINT_ROOT_CONCAT(IQTREE_CONCATENATED.out.tree)
    }

    // --- Step 7: Barcode alignment visualization (always with OUTGROUP as reference) ---
    VISUALIZE_ALIGNMENT(
        CONCATENATE_ALIGNMENTS.out.fasta,
        CONCATENATE_ALIGNMENTS.out.partition
    )
}
