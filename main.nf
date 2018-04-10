#!/usr/bin/env nextflow

/* 
 * Main embeded analysis pipeline script
 *
 * @authors
 * Evan Floden <evanfloden@gmail.com> 
 */

params.name             = "emebed-analysis-nf"
params.ref              = "$baseDir/tutorial/refs/seatoxin.ref"
params.seqs             = "$baseDir/tutorial/seqs/seatoxin.*.*.fa"
params.output           = "$baseDir/results/"
params.default_methods  = "CLUSTALO"
params.std_methods      = "CLUSTALO_STD"
params.dpa_methods      = "CLUSTALO_DPA"
params.tree_methods     = "CLUSTALO"
params.buckets          = '250,1000'

// create dpa alignments [BOOL]
params.dpa_align = true

// create standard alignments [BOOL]
params.std_align = true

// create default alignments [BOOL]
params.default_align = true

log.info "e m b e d e d  -  a n a l y s i s  ~  version 0.2"
log.info "====================================="
log.info "name                            : ${params.name}"
log.info "sequences (FA)                  : ${params.seqs}"
log.info "reference alignment (ALN)       : ${params.ref}"
log.info "output (DIRECTORY)              : ${params.output}"
log.info "aligners [CLUSTALO|MAFFT]       : ${params.default_methods} | ${params.std_methods} | ${params.dpa_methods}"
log.info "tree methods                    : ${params.tree_methods}"
log.info "bucket sizes                    : ${params.buckets}"
log.info "\n"

/**************************
 * 
 * S E T U P   I N P U T S   A N D   P A R A M E T E R S
 *
 */

/*
 * Create a channel for input alignment files & the seeds id files
 */

Channel
    .fromPath( params.ref )
    .ifEmpty { error "Cannot find any input reference files matching: ${params.ref}" }
    .map { file -> tuple( file.baseName, file ) }
    .into { refs1; refs2 }

Channel
    .fromPath( params.seqs )
    .ifEmpty { error "Cannot find any input files matching: ${params.seqs}" }
    .map { file -> tuple( file.baseName, file ) }
    .map { filename,  file -> def (dataset, size, replicate) = filename.tokenize('.'); [dataset, size, replicate, file] } 
    .view()
    .into { seqs1; seqs2; seqs3; seqs4}
 
refs1
    .cross(seqs1)
                    //   datset    size       replicate    refFile    combinedSeqFile
    .map { item -> [item[0][0], item[1][1], item[1][2], item[0][1], item[1][1]] }
    .set { seqsAndRefs }

tree_methods = params.tree_methods

/*
 *
 **************************/


process guide_trees {
   tag "${id}.${tree_method}.${size}.${rep}"
   publishDir "${params.output}/guide_trees", mode: 'copy', overwrite: true

   input:
     set val(id),   \
         val(size), \
         val(rep),  \
         file(seqs) \
         from seqs2

     each tree_method from tree_methods.tokenize(',')

   output:
     set val(id), 
         val(size),
         val(rep),
         val(tree_method),
         file("${id}.${tree_method}.${size}.${rep}.dnd") \
         into treesGenerated

   when:
     params.std_align || params.dpa_align

   script:
     template "tree/generate_tree_${tree_method}.sh"
}



treesGenerated
    //                [id  size,  rep],   id,  size,    rep, treeMethod, treeFile
    .map { it -> [ [it[0],it[1],it[2]], it[0], it[1], it[2], it[3], it[4]] }
    .set { treesMapped }

seqs3
    .map { it -> [ [it[0],it[1],it[2]], it[0], it[1], it[2], it[3] ] }
    .set { sequencesMapped }

treesMapped
    .combine (sequencesMapped, by:0)
    .map { it -> [it[1], it[2], it[3], it[4], it[5], it[9] ]}
    .into {sequenceSetsWithTreesForDPA;sequenceSetsWithTreesForSTD}



process std_alignment {

    tag "${id}.${size}.${rep}.${align_method}.STD.NA.${tree_method}"
    publishDir "${params.output}/alignments", mode: 'copy', overwrite: true

    input: 
      set val(id),         \
          val(size),       \
          val(rep),        \
          val(tree_method),\
          file(guide_tree),\
          file(seqs)       \
          from sequenceSetsWithTreesForSTD

      each align_method from params.std_methods.tokenize(',')

    when:
      params.std_align

    output:
      set val(id), 
          val(size),
          val(rep),
          val(align_method),
          val(tree_method), 
          val("std_align"), 
          val("NA"), 
          file("${id}.${size}.${rep}.std.NA.${align_method}.with.${tree_method}.tree.aln") \
      into std_alignments

     script:
       template "std_align/std_align_${align_method}.sh"
}


process dpa_alignment {

    tag "${id}.${size}.${rep}.${align_method}.DPA.${bucket_size}.${tree_method}"
    publishDir "${params.output}/alignments", mode: 'copy', overwrite: true

    input:
      set val(id),         \
          val(size),       \
          val(rep),        \
          val(tree_method),\
          file(guide_tree),\
          file(seqs)       \
          from sequenceSetsWithTreesForDPA

      each bucket_size from params.buckets.tokenize(',')

      each align_method from params.dpa_methods.tokenize(',')

    output:
      set val(id),
          val(size),
          val(rep),
          val("${align_method}"), \
          val(tree_method), \
          val("dpa_align"), \
          val(bucket_size), \
          file("${id}.${size}.${rep}.dpa.${bucket_size}.${align_method}.with.${tree_method}.tree.aln") \
          into dpa_alignments

    when:
      params.dpa_align

    script:
       template "dpa_align/dpa_align_${align_method}.sh"
}



process default_alignment {

    tag "${id}.${size}.${rep}.${align_method}.default.NA.default" 
   publishDir "${params.output}/alignments", mode: 'copy', overwrite: true

    input:
      set val(id),   \
          val(size), \
          val(rep),  \
          file(seqs) \
          from seqs4

      each align_method from params.default_methods.tokenize(',')

    when:
      params.default_align

    output:
      set val(id), 
          val(size),
          val(rep),
          val("${align_method}"), \
          val("default"), \
          val("default_align"), \
          val("NA"), \
          file("${id}.${size}.${rep}.default.NA.${align_method}.with.default.tree.aln") \
          into default_alignments

     script:
       template "default_align/default_align_${align_method}.sh"
}



std_alignments
  .mix ( dpa_alignments )
  .mix (default_alignments )
  .set { all_alignments }

all_alignments
    .combine(refs2, by:0) 
    .set { toEvaluate } 

// [ val(id), 
//   val(size), 
//   val(rep),
//   val(align_method),
//   val(tree_method), 
//   val(align_type),
//   val(bucket_size),
//   file (alignment), 
//   file(reference) ]


process evaluate {

    tag "${id}.${size}.${rep}.${align_method}.${align_type}.${bucket_size}.${tree_method}"

    input:
      set val(id), 
          val(size),              \
          val(rep),               \
          val(align_method),      \
          val(tree_method),       \
          val(align_type),        \
          val(bucket_size),       \
          file (test_alignment),  \
          file(ref_alignment)     \
          from toEvaluate

    output:
      set val(id),
          val(size),
          val(rep), 
          val(align_method),
          val(align_type),
          val(bucket_size),
          val(tree_method),
          file("score.sp.tsv") \
          into spScores

      set val(id), 
          val(size),
          val(rep),
          val(align_method),
          val(align_type),
          val(bucket_size),
          val(tree_method),
          file("score.tc.tsv") \
          into tcScores

      set val(id), 
          val(size),
          val(rep),
          val(align_method),
          val(align_type),
          val(bucket_size),
          val(tree_method),
          file("score.col.tsv") \
          into colScores

     script:
     """
       t_coffee -other_pg aln_compare \
                -al1 ${ref_alignment} \
                -al2 ${test_alignment} \
                -compare_mode sp \
            | grep -v "seq1" |grep -v '*' | awk '{ print \$4}' ORS="\t" \
            >> "score.sp.tsv"

       t_coffee -other_pg aln_compare \
                -al1 ${ref_alignment} \
                -al2 ${test_alignment} \
                -compare_mode tc \
            | grep -v "seq1" |grep -v '*' | awk '{ print \$4}' ORS="\t" \
            >> "score.tc.tsv"

       t_coffee -other_pg aln_compare \
                -al1 ${ref_alignment} \
                -al2 ${test_alignment} \
                -compare_mode column \
            | grep -v "seq1" |grep -v '*' | awk '{ print \$4}' ORS="\t" \
            >> "score.col.tsv"
    """
}


spScores
    .collectFile(name:"spScores.${workflow.runName}.csv", newLine:true, storeDir: "$params.output/scores" ) { 
        it[0]+"\t"+it[1]+"\t"+it[2]+"\t"+it[3]+"\t"+it[4]+"\t"+it[5]+"\t"+it[6]+"\t"+it[7].text }

tcScores
    .collectFile(name:"tcScores.${workflow.runName}.csv", newLine:true, storeDir: "$params.output/scores" ) {
        it[0]+"\t"+it[1]+"\t"+it[2]+"\t"+it[3]+"\t"+it[4]+"\t"+it[5]+"\t"+it[6]+"\t"+it[7].text }

colScores
    .collectFile(name:"colScores.${workflow.runName}.csv", newLine:true, storeDir: "$params.output/scores" ) {
        it[0]+"\t"+it[1]+"\t"+it[2]+"\t"+it[3]+"\t"+it[4]+"\t"+it[5]+"\t"+it[6]+"\t"+it[7].text }

/*
 *
 **************************/
