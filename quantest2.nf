
// Dir
params.maindir = "/users/cn/sjin/projects/homoplasy/slaveTree_results"  // Where input & output

// Input alignments
params.msa = "${params.maindir}/results_*_MAFFT-GINSI/alignments/*.aln"

// Input fasta sequences
params.seqs = "/users/cn/egarriga/datasets/homfam/combinedSeqs/*.fa"

// Input references: aux & ss
params.ref_aux = "/users/cn/sjin/projects/homoplasy/ss/informative3/*.aux"
params.ref_ss = "/users/cn/sjin/projects/homoplasy/ss/informative3/*.ss"

// // Alignment bucket & method
// params.bucket = "50,100,200,500,1000"
// params.aligner = "CLUSTALO"   // MAFFT-FFTNS1,MAFFT-SPARSECORE,MAFFT-GINSI

// which guide tree to retrieve N informative sequences 
params.regtrim_tree = "/users/cn/sjin/projects/homoplasy/trees/*.mafftdnd.dnd"
params.n_quantest2 = "1000"


// ----------------
// Print input info
// ----------------

log.info """\
         Running QuanTest2 - modified version"
         ======================================="
         Main directory                                     : ${params.maindir} 
         MSA                                                : ${params.msa}        
         Sequences                                          : ${params.seqs}
         References aux                                     : ${params.ref_aux}
         References ss                                      : ${params.ref_ss}
         Guide tree to retrieve informative sequences       : ${params.regtrim_tree}
         N seq --> QuanTest2                                : ${params.n_quantest2}
         """
         .stripIndent()


// ---------------------------------
// Orgaize input files into Channels
// ---------------------------------

Channel
  .fromPath( params.msa )
  .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[5], 
                   item.baseName.tokenize('.')[2], item.baseName.tokenize('.')[3], 
                   item.baseName, item] }   // [ fam, tree_method, bucket, aligner, id, msa]
  .set { ch_msa }

Channel
  .fromPath( params.seqs )
  .map { item -> [ item.baseName, item] }  
  .set { ch_seqs }

Channel
  .fromPath( params.ref_aux )
  .map { item -> [ item.baseName, item] }
  .set { ch_refaux }

Channel
  .fromPath( params.ref_ss )
  .map { item -> [ item.baseName, item] }
  .set { ch_refss }

Channel
  .fromPath( params.regtrim_tree )
  .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }   //  [fam, regtrimtree_method, regtrimtree]
  .set { ch_regtrimtree }

ch_msa
  .combine( ch_seqs, by:0 )
  .combine( ch_refaux, by:0 )
  .combine( ch_refss, by:0 )
  .combine( ch_regtrimtree, by:0 )
  .set { toQuanTest2 }
//   .println()


// -------------
// Run QuanTest2
// -------------

process quantest2 {

    // cache false 

    tag "${id}.informative${n_quantest2}.with.${regtrim_treemethod}.tree*"
    publishDir "${params.maindir}/results_${bucket_size}_${aligner_method}/quantest2", mode: 'copy', overwrite: true

    input:
      set val(fam), \
          val(tree_method), \
          val(bucket_size), \
          val(aligner_method), \
          val(id), \
          file(msa), \
          file(seqs), \
          file(refaux), \
          file(refss), \
          val(regtrim_treemethod), \
          file(regtrim_tree) \
          from toQuanTest2        
      each n_quantest2 from params.n_quantest2.tokenize(',')
          
    output:
      set val(id), \
          val(fam), \
          val(tree_method), \
          val(bucket_size), \
          val(aligner_method), \
          val(regtrim_treemethod), \
          val(n_quantest2), \
          file("${id}.informative${n_quantest2}.with.${regtrim_treemethod}.tree.quantest2*") \
          into predictions

    // when:
    //   if(new File("${params.output}/n${n_quantest2}/${fam}/${id}.informative${n_quantest2}.with.${regtrim_tree}.tree.quantest2").exists()){false}else{true}

    script:
      """
      python3 $baseDir/bin/get_refANDinformative_seqs.py \
      --msa ${msa} --seq ${seqs} --names ${refaux} --tree ${regtrim_tree} --n ${n_quantest2} \
      >> ${id}.informative${n_quantest2}.with.${regtrim_treemethod}.tree.aln




      quantest2 ${id}.informative${n_quantest2}.with.${regtrim_treemethod}.tree.aln ${refss} 

      awk 'NR==4{print}' ${id}.informative${n_quantest2}.with.${regtrim_treemethod}.tree.quantest2_log | cut -f2 \
      > ${id}.informative${n_quantest2}.with.${regtrim_treemethod}.tree.quantest2
      
      """

}


workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}" 
}
