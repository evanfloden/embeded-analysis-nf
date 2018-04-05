#!/usr/bin/env nextflow

/* 
 * Main embeded analysis pipeline script
 *
 * @authors
 * Evan Floden <evanfloden@gmail.com> 
 */

params.name             = "emebed-analysis-nf"
params.ref              = "$baseDir/tutorial/seatoxin.ref"
params.seqs             = "$baseDir/tutorial/seatoxin.fa"
params.output           = "$baseDir/results/"
params.reps             = 10
params.align_method     = "CLUSTALO"
params.tree_method      = "CLUSTALO"
params.buckets          = '250,1000'

log.info "e m b e d e d  -  a n a l y s i s  ~  version 0.1"
log.info "====================================="
log.info "name                            : ${params.name}"
log.info "test sequences (FA)             : ${params.seqs}"
log.info "reference alignment (ALN)       : ${params.ref}"
log.info "output (DIRECTORY)              : ${params.output}"
log.info "aligners [CLUSTALO|MAFFT]       : ${params.align_method}"
log.info "tree methods                    : ${params.tree_method}"
log.info "bucket sizes                    : ${params.buckets}"
log.info "replicates                      : ${params.reps}"
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
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.ref}" }
    .map { file -> tuple( file.baseName, file ) }
    .into { refs1; refs2 }

Channel
    .fromPath( params.seqs )
    .ifEmpty { error "Cannot find any input seeds files matching: ${params.seeds}" }
    .map { file -> tuple( file.baseName, file ) }
    .into { seqs1; seqs2 }

refs1
    .cross(seqs1)
    .map { item -> [item[0][0], item[0][1], item[1][1]] }
    .set { seqsAndRefs }

/*
 *
 **************************/


/**************************
 *
 *   C A L C U L A T E   S E Q  N U M B E R S  F O R  
 *   E A C H  D A T A S E T / A L I G N M E N T
 *
 */

process calculate_seqs {

    input:
        set val(datasetID), val(seqs) from seqs2

    output:
        set val(datasetID), val(number_seqs_per_alignment) into seqNumbers

    exec:
    number_seqs_per_alignment = [0,5,10,25,50,100,200,400,600,800,1000,1250,\
                                1500,1750,2000,2500,3000,4000,5000,7500,\
                                10000,15000,20000,50000,100000]

}


def transform = { 
          def result = []
          def name = it[0]
          it[1].each { result << [name, it] }
          return result
 }


seqNumbers
    .flatMap(transform)
    .set {seqNumbersFlat}

seqsAndRefs
    .cross(seqNumbersFlat)
    .map { item -> [ item[1][0], item[1][1], item[0][1], item[0][2] ] }
    .set { seqsAndRefsAndNumbers }

/**************************
 *
 *   G E N E R A T E   S E Q U E N C E   S E T S  T O   A L I G N
 *
 */

process generate_sequence_sets {

    publishDir "${params.output}/sequence_sets/${datasetID}/", mode: 'copy';
    tag "${datasetID} - ${size} - ${rep}"

    input:
        set val(datasetID), val(size), file(references), file(sequences) from seqsAndRefsAndNumbers
        each rep from (1..params.reps)

    output:
        set val(datasetID), val(size), val(rep), file ("${datasetID}.${size}.${rep}.fa") into sequenceSets, sequenceSets2

    script:
    """
    # FORMAT A FASTA FILE CONTAINING ALL REF SEQUENCES
    t_coffee -other_pg seq_reformat -in ${references} -output fasta_seq -out refs.tmp.fa

    # FORMAT A FASTA FILE CONTAINING ALL OTHER SEQUENCES
    esl-reformat fasta ${sequences} > seqs.tmp.fa

    # CREATE FILE OF SEQ IDS (HEADERS), SORTED RANDOMLY
    grep '^>' seqs.tmp.fa | sort -R > headers.txt

    # SELECT FIRST ${size} HEADERS
    head -n ${size} headers.txt > headers.${size}.txt
    sed -i 's/^.\\{1\\}//g' headers.${size}.txt 

    # SELECTS SEQS IN headers.${size}.txt
    faSomeRecords ${sequences} headers.${size}.txt selected_seqs.fa

    # COMBINE SELECTED SEQS AND REFS
    cat selected_seqs.fa >> refs.tmp.fa 

    # FORMAT AND SHUFFLE
    t_coffee -other_pg seq_reformat \
             -in refs.tmp.fa \
             -output fasta_seq \
             -out ${datasetID}.${size}.${rep}.fa \
             -action +reorder random
    """
}

process guide_trees {
   tag "${datasetID} - ${tree_method} - ${size} - ${rep}"
   publishDir "${params.output}/guide_trees", mode: 'copy', overwrite: true

   input:
     set val(datasetID),  
         val(size), 
         val(rep), 
         file(sequenceSet) from sequenceSets2

     each tree_method from params.tree_method.tokenize(',') 

   output:
     set val(datasetID), \
         val(tree_method), 
         val(size), \
         val(rep), \
         file("${datasetID}.${tree_method}.${size}.${rep}.dnd") \
         into treesGenerated

   script:
     template "trees/generate_tree_${tree_method}.sh"
}


treesGenerated
    .map { it -> [ [it[0],it[2],it[3]], it[0], it[1], it[2], it[3], it[4]] }
    .set { treesMapped }

sequenceSets
    .map { it -> [ [it[0],it[1],it[2]], it[0], it[1], it[2], it[3] ] }
    .set { sequencesMapped }

treesMapped
    .combine (sequencesMapped, by:0)
    .map { it -> [it[1],it[2], it[5], it[3], it[4], it[9] ]}
    .set {sequenceSetsWithTrees}

process generate_dpa_alignments {

    tag "${datasetID} - ${align_method} - ${tree_method} - DPA - ${bucket_size} - ${size} - ${rep}"

    publishDir "${params.output}/alignments/${datasetID}/", mode: 'copy';
    
    input:
        set val(datasetID), val(tree_method), file(tree), val(size), val(rep), file (seqs2align) from sequenceSetsWithTrees     
        
        each bucket_size from params.buckets.tokenize(',')
 
        each align_method from params.align_method.tokenize(',') 

   output:
        set val(datasetID), val(align_method), val(tree_method), val(bucket_size), val(size), val(rep), file ("${datasetID}.${align_method}.${tree_method}.${bucket_size}.${size}.${rep}.aln") into completeAlignments
    
    script:
    template "align/generate_alignments_dpa_${align_method}.sh"
}   


completeAlignments
    .combine(refs2, by:0) 
    .set { toEvaluate } 

// [ val(datasetID), val(align_method), val(tree_method), val(size), val(rep), file (alignment), file(ref) ]


process evaluate {

    tag "${id} - ${tree_method} - ${align_method} - ${size} - ${rep} - ${bucket_size}"

    input:
      set val(id), \
          val(align_method), \
          val(tree_method), \
          val(bucket_size),\
          val(size), \
          val(rep), \
          file(test_alignment), \
          file(ref_alignment) \
          from toEvaluate

    output:
      set val(id), val(align_method), \
          val(tree_method), val(bucket_size), \
          val(size), val(rep), file("score.sp.tsv") \
          into spScores

      set val(id), val(align_method), \
          val(tree_method), val(bucket_size), \
          val(size), val(rep), file("score.tc.tsv") \
          into tcScores

       set val(id), val(align_method), \
          val(tree_method), val(bucket_size), \
          val(size), val(rep),file("score.col.tsv") \
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
        it[0]+"\t"+it[1]+"\t"+it[2]+"\t"+it[3]+"\t"+it[4]+"\t"+it[5]+"\t"+it[6].text }

tcScores
    .collectFile(name:"tcScores.${workflow.runName}.csv", newLine:true, storeDir: "$params.output/scores" ) {
        it[0]+"\t"+it[1]+"\t"+it[2]+"\t"+it[3]+"\t"+it[4]+"\t"+it[5]+"\t"+it[6].text }

colScores
    .collectFile(name:"colScores.${workflow.runName}.csv", newLine:true, storeDir: "$params.output/scores" ) {
        it[0]+"\t"+it[1]+"\t"+it[2]+"\t"+it[3]+"\t"+it[4]+"\t"+it[5]+"\t"+it[6].text }

/*
 *
 **************************/
