/* 
 * Main embeded analysis pipeline script
 *
 * @authors
 * Evan Floden <evanfloden@gmail.com> 
 * Edgar
 */

params.name             = "emebed-analysis-nf"
params.ref              = "$baseDir/tutorial/seatoxin.ref"
params.seqs             = "$baseDir/tutorial/seatoxin.fa"
params.output           = "$baseDir/results/"
params.alignments       = 50
params.aligner          = "CLUSTALO"

log.info "e m b e d e d  -  a n a l y s i s  ~  version 0.1"
log.info "====================================="
log.info "name                            : ${params.name}"
log.info "test sequences (FA)             : ${params.seqs}"
log.info "reference alignment (ALN)       : ${params.ref}"
log.info "output (DIRECTORY)              : ${params.output}"
log.info "number of alignments            : ${params.alignments}"
log.info "aligner [CLUSTALO|MAFFT|UPP]    : ${params.aligner}"
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
    .view()
    .set { seqsAndRefs }

/*
 * Select aligner to generate MSA:
 *
 * 'CLUSTALO' | 'MAFFT' |  'UPP'
 *
 */
    aligner = params.aligner

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
    alignments = params.alignments as int
    min = 1
    seqFileText = seqs.text
    def max = 0
    seqFileText.eachLine { String line ->
        if(line.startsWith('>')) {
            max++ 
        }   
    }   
    println "max = $max"
    
    // Define linear function for scale
    number_seqs_per_alignment = []
    
    // Define the number of sequences in each alignment
    number_seqs_per_alignment.add(0)
    number_seqs_per_alignment.add(min)
    for (x = 1; x < alignments-1; x++) {
      k = (max-min)/(alignments-1)*x
      l = k.intValue()
      number_seqs_per_alignment.add(l)
    } 
    number_seqs_per_alignment.add(max)

}

seqNumbers
    .flatMap { name = it[0]; return it[1] }
    .map { item -> [name, item] }
    .set {seqNumbersFlat}

seqsAndRefs
    .cross(seqNumbersFlat)
    .map { item -> [ item[1][0], item[1][1], item[0][1], item[0][2], item[1][1] ] }
    .set { seqsAndRefsAndNumbers }

/**************************
 *
 *   G E N E R A T E   S E Q U E N C E   S E T S  T O   A L I G N
 *
 */

process generate_sequence_sets {
    container "cbcrg/tcoffee"  
    tag "${datasetID}_${size}_${rep}"

    input:
        set val(datasetID), val(size), file(alnFile), file(seqs) from seqsAndRefsAndNumbers
        each rep from (0..9)

    output:
        set val(datasetID), val(size), val(rep), file ("${size}.${rep}.fa") into sequenceSets

    script:
    """
    grep '^>' ${seqs} | sort -R > headers.txt
    head -n ${size} headers.txt > headers.${size}.txt
    sed 's/^.\\{1\\}//g' headers.${size}.txt > headers.${size}.x.txt
    faSomeRecords ${seqs} headers.${size}.x.txt shuffled_seqs.fa
    esl-reformat --informat afa fasta ${alnFile} > ref_seqs.fa
    cat ref_seqs.fa > ${size}.${rep}.fa
    cat shuffled_seqs.fa >> ${size}.${rep}.fa  
    """
}


process generate_alignments {
    container "cbcrg/benchfam_large_scale"
    tag "${datasetID}_${size}_${rep}"
    publishDir "${params.output}/${aligner}/${datasetID}/alignments/${size}/", mode: 'copy';
    
    input:
        set val(datasetID), val(size), val(rep), file ("${size}.${rep}.fa") from sequenceSets       
 
    output:
        set val(datasetID), val(size), val(rep), file ("${size}.${rep}.afa") into completeAlignments
    
    script:
    template "generate_alignments_${aligner}.sh"
}   




refs2
    .cross(completeAlignments) 
    .map { item -> [ [item[0][0], item[1][1]], item[1][3], item[0][1]] }
    .groupTuple()
    .map { item -> [ item[0][0], item[0][1], item[1], item[2][0]] }
    .set{ alignmentsGrouped }



process evaluate_alignments {

    container "cbcrg/tcoffee"
    tag "${dataset} ${size}"

    input:
    set val(dataset), val(size), file(alignments), file(refAln) from alignmentsGrouped

    output:
    set val(dataset), val(size), file("${aligner}.${dataset}.${size}.sp") into spScores
    set val(dataset), val(size), file("${aligner}.${dataset}.${size}.tc") into tcScores
    set val(dataset), val(size), file("${aligner}.${dataset}.${size}.col") into colScores

    script:
    """
    touch ${aligner}.${dataset}.${size}.sp
    for i in {0..9};
    do 
    t_coffee -other_pg aln_compare \
             -al1 ${refAln} \
             -al2 ${size}.\$i.afa \
            -compare_mode sp \
            | grep -v "seq1" |grep -v '*' | awk '{ print \$4}' ORS="\t" \
            >> "${aligner}.${dataset}.${size}.sp"
    done

    touch ${aligner}.${dataset}.${size}.tc
    for i in {0..9};
    do
    t_coffee -other_pg aln_compare \
             -al1 ${refAln} \
             -al2 ${size}.\$i.afa \
            -compare_mode tc \
            | grep -v "seq1" |grep -v '*' | awk '{ print \$4}' ORS="\t" \
            >> "${aligner}.${dataset}.${size}.tc"
    done

    touch ${aligner}.${dataset}.${size}.col
    for i in {0..9};
    do
    t_coffee -other_pg aln_compare \
             -al1 ${refAln} \
             -al2 ${size}.\$i.afa \
            -compare_mode column \
            | grep -v "seq1" |grep -v '*' | awk '{ print \$4}' ORS="\t" \
            >> "${aligner}.${dataset}.${size}.col"
    done

   """
}

spScores
    .collectFile(name:"spScores.csv", sort:{ it[0] }, newLine:true, storeDir: params.output ) { 
        it[0]+"\t"+it[1]+"\t"+it[2].text }

tcScores
    .collectFile(name:"tcScores.csv", sort:{ it[0] }, newLine:true, storeDir: params.output ) {
        it[0]+"\t"+it[1]+"\t"+it[2].text }

colScores
    .collectFile(name:"colScores.csv", sort:{ it[0] }, newLine:true, storeDir: params.output ) {
        it[0]+"\t"+it[1]+"\t"+it[2].text }

/*
 *
 **************************/
