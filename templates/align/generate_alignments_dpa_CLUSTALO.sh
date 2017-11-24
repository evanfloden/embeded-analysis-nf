t_coffee -dpa -dpa_method clustalo_msa \
         -dpa_tree ${tree} \
         -seq ${seqs2align} \
         -dpa_nseq ${bucket_size} \
         -outfile ${datasetID}.${align_method}.${tree_method}.${bucket_size}.${size}.${rep}.aln
