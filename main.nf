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
params.alignments       = 3
params.min              = 1
params.max              = 50
params.aligner          = "CLUSTALO"

log.info "e m b e d e d  -  a n a l y s i s  ~  version 0.1"
log.info "====================================="
log.info "name                            : ${params.name}"
log.info "test sequences (FA)             : ${params.seqs}"
log.info "reference alignment (ALN)       : ${params.ref}"
log.info "output (DIRECTORY)              : ${params.output}"
log.info "number of alignments            : ${params.alignments}"
log.info "minimum number of seqs          : ${params.min}"
log.info "minimum number of seqs          : ${params.max}"
log.info "aligner [CLUSTALO|MAFFT|UPP]    : ${params.aligner}"
log.info "\n"


/**************************
 * 
 * S E T U P   I N P U T S   A N D   P A R A M E T E R S
 *
 */
    alignments             = params.alignments as int
    min                    = params.min as int
    max                    = params.max as int
  
    // Define linear function for scale
    number_seqs_per_alignment = []

    // Define the number of sequences in each alignment
    number_seqs_per_alignment.add(0)
    number_seqs_per_alignment.add(min)
    for (x = 1; x < alignments-1; x++) {
      k = (max-min)/(alignments-1)*x 
      BigInteger l = BigInteger.valueOf(k.intValue());
      number_seqs_per_alignment.add(l)
    }
    number_seqs_per_alignment.add(max)

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
    .into { seqs }

refs1
    .cross(seqs)
    .map { item -> [item[0][0], item[0][1], item[1][1]] }
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
 *   G E N E R A T E   A L I G N M E N T S
 *
 */

process generate_alignments {
 
    tag "${datasetID}_${x}_${rep}"
    publishDir "${params.output}/${datasetID}/alignments/${x}/${rep}", mode: 'copy';

    input:
        set val(datasetID), file(alnFile), file(seqs) from seqsAndRefs
        each x from number_seqs_per_alignment
        each rep from (0..9)

    output:
        set val(datasetID), val(x), val(rep), file ("${x}.${rep}.fa") into completeAlignments

    script:
    """

    grep '^>' ${seqs} | sort -R > headers.txt

    head -n ${x} headers.txt > headers.${x}.txt

    sed 's/^.\\{1\\}//g' headers.${x}.txt > headers.${x}.x.txt

    faSomeRecords ${seqs} headers.${x}.x.txt shuffled_seqs.fa

    esl-reformat --informat afa fasta ${alnFile} > ref_seqs.fa

    cat ref_seqs.fa > seqs2align.${x}.fa
    cat shuffled_seqs.fa >> seqs2align.${x}.fa
   
    clustalo --infile=seqs2align.${x}.fa --outfmt=fa --force -o ${x}.${rep}.fa  

    """
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
    publishDir "${params.output}/${aligner}/${datasetID}/evaluation/", mode: 'copy';

    input:
    set val(dataset), val(size), file(alignments), file(refAln) from alignmentsGrouped

    output:
    set val(dataset), val(size), file("${aligner}.${dataset}.${size}.sp") into spScores
    set val(dataset), val(size), file("${aligner}.${dataset}.${size}.sp") into tcScores
    set val(dataset), val(size), file("${aligner}.${dataset}.${size}.sp") into colScores

    script:
    """
    touch ${aligner}.${dataset}.${size}.sp
    for i in {0..9};
    do 
    t_coffee -other_pg aln_compare \
             -al1 ${refAln} \
             -al2 ${size}.\$i.fa \
            -compare_mode sp \
            | grep -v "seq1" |grep -v '*' | awk '{ print \$4}' ORS="\t" \
            >> "${aligner}.${dataset}.${size}.sp"
    done

    touch ${aligner}.${dataset}.${size}.tc
    for i in {0..9};
    do
    t_coffee -other_pg aln_compare \
             -al1 ${refAln} \
             -al2 ${size}.\$i.fa \
            -compare_mode tc \
            | grep -v "seq1" |grep -v '*' | awk '{ print \$4}' ORS="\t" \
            >> "${aligner}.${dataset}.${size}.tc"
    done

    touch ${aligner}.${dataset}.${size}.col
    for i in {0..9};
    do
    t_coffee -other_pg aln_compare \
             -al1 ${refAln} \
             -al2 ${size}.\$i.fa \
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
