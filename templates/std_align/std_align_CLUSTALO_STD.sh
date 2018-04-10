clustalo --infile=${seqs} \
         --guidetree-in=${guide_tree} \
         --outfmt=fa \
         -o ${id}.${size}.${rep}.std.NA.${align_method}.with.${tree_method}.tree.aln
