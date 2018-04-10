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

log.info "e m b e d e d  -  a n a l y s i s  ~  version 0.1"
log.info "====================================="
log.info "name                            : ${params.name}"
log.info "test sequences (FA)             : ${params.seqs}"
log.info "reference alignment (ALN)       : ${params.ref}"
log.info "output (DIRECTORY)              : ${params.output}"
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

    publishDir "${params.output}/sequence_sets/", mode: 'copy';
    tag "${datasetID} - ${size} - ${rep}"

    input:
        set val(datasetID), val(size), file(references), file(sequences) from seqsAndRefsAndNumbers
        each rep from (1..params.reps)

    output:
        set val(datasetID), val(size), val(rep), file ("${datasetID}.${size}.${rep}.fa") into sequenceSets

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
  
     # Make null if number seqs less than size 
     x=\$(grep '>' ${datasetID}.${size}.${rep}.fa | wc -l) 
     if (( \$x < ${size} )); then
       echo 'null' > ${datasetID}.${size}.${rep}.fa
     fi 
 
    """
}

