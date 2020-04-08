#!/usr/bin/env nextflow



// which alignment methods to run
params.bucket = [100,1000]
params.aligner = ["CLUSTALO":"CO", "MAFFT-FFTNS1":"FFTNS1"] //CLUSTALO,MAFFT-FFTNS1,MAFFT-SPARSECORE,UPP,MAFFT-GINSI"
params.tree = ["codnd","FAMSA","mafftdnd","dpparttreednd1","fastaparttreednd","fftns1dnd","parttreednd0","parttreednd2","CLUSTALO-RANDOM"]
//params.tree = ["codnd","FAMSA"]

// which guide tree to retrieve informative sequences 
params.regtrim_tree = "mafftdnd"
params.n_quantest2 = [1000,500,300,250,200,150,100,90,80,70,60,50,40,30]  // Size of the partial alignment given to QuanTest2
//params.n_quantest2 = [100, 1000]

// which families to use
//params.fams = [ "sdr", "Acetyltransf" ]
//params.fams = ["sdr","aadh"]
//params.fams=[ "ltn", "aadh" ]
params.fams = ["sdr","Acetyltransf","rrm","aat","adh","p450","rhv","blmb","PDZ","Rhodanese","hla","aldosered","ghf13","hom","biotin_lipoyl","tRNA-synt_2b","myb_DNA-binding","gluts","blm","egf","gpdh","lyase_1","int","subt","ldh","HLH","LIM","cyclo","proteasome","icd","msb","OTCace","HMG_box","flav","uce","peroxidase","sodfe","ghf1","cys","ace","glob","tim","hr","hormone_rec","hpr","oxidored_q6","asp","cytb","serpin","annexin","aadh","phc","ghf5","Ald_Xan_dh_2","mofe","Sulfotransfer","kunitz","GEL","tms","DMRL_synthase","KAS","sodcu","tgfb","ghf10","rub","mmp","cah","DEATH","cryst","kringle","az","il8","ltn"]

//

// input sequence directory
params.seq_dir = "/users/cn/egarriga/datasets/homfam/combinedSeqs"  // *.fa

// input reference sequences (3 per dataset) directory
params.ref_dir = "/users/cn/sjin/projects/homoplasy/ss/informative3"  // *.aux

// input reference secondary structure (3 per dataset) directory
params.ss_dir = "/users/cn/sjin/projects/homoplasy/ss/informative3"  // *.ss

// input guide tree (for T-coffee +regtrim) directory
params.tree_dir = "/users/cn/sjin/projects/homoplasy/trees"   // *.dnd

// input msa directory
params.msa_dir = "/users/cn/sjin/projects/homoplasy/nf_homoplasty/"   // subfolder 'results_fullTree_ + .... '

// output directory
params.output = "/users/cn/sjin/projects/homoplasy/test2"

//

log.info """\
         Running QuanTest2 - modified version"
         ======================================="
         Dataset                                            : ${params.fams}
         Alignment generated with:
            Bucket size                                     : ${params.bucket}
            Aligner                                         : ${params.aligner}
            Guide tree                                      : ${params.tree}
         Guide tree to retrieve informative sequences       : ${params.regtrim_tree}
         N seq --> QuanTest2                                : ${params.n_quantest2}
         Output directory (DIRECTORY)                       : ${params.output}
         """
         .stripIndent()





// organize input files 

msas_l=[]; seqs_l=[]; names_l=[]; trees_l=[]; ss_l=[]
for ( fam in params.fams ) 
{
  for ( bucket in params.bucket ) 
  {
    for ( aligner in params.aligner.keySet().collect() ) 
    {
        f=params.aligner[aligner]
        for ( tree in params.tree ) 
        {
          msa=params.msa_dir+"/results_fullTree_"+f+"_bucket"+bucket+"/alignments/"+fam+".dpa_"+bucket+"."+aligner+".with."+tree+".tree.aln"
          s=params.seq_dir + "/" + fam + ".fa"
          n=params.ref_dir + "/" + fam + ".aux"
          t=params.tree_dir + "/" + fam + "." + params.regtrim_tree + ".dnd"
          ss=params.ss_dir + "/" + fam + ".ss"
          msas_l.add(msa); seqs_l.add(s); names_l.add(n); trees_l.add(t); ss_l.add(ss)
        }
    }
  }
}


Channel
  .fromPath(msas_l)
  .map { item -> [ item.baseName ]}
  .set { ids }

Channel
  .fromPath(seqs_l)
  .map { item -> [ item.baseName ]}
  .set { fams }

Channel
  .fromPath(msas_l)
  .map { item -> [ item ]}
  .set { msas }

Channel
  .fromPath(seqs_l)
  .map { item -> [ item ]}
  .set { seqs }

Channel
  .fromPath(names_l)
  .map { item -> [ item ]}
  .set { refnames }

Channel
  .fromPath(trees_l)
  .map { item -> [ item ]}
  .set { trees }

Channel
  .fromPath(ss_l)
  .map { item -> [ item ]}
  .set { ss }

ids
  .merge( fams, msas, seqs, refnames, trees, ss )
  .set { toQuantest2 }   // Channel [ id, fam, msa, seq, refname, tree, ss ]




process nf_quantest2 {

    tag "${id}.informative${n_quantest2}.with.${regtrim_tree}.tree"
    publishDir "${params.output}/n${n_quantest2}/${fam}", mode: 'copy', overwrite: true

    input:
      set val(id), \
          val(fam), \
          file(msa), \
          file(seqs), \
          file(refnames), \
          file(tree), \
          file(ss) \
          from toQuantest2        
      val(regtrim_tree) from params.regtrim_tree
      each n_quantest2 from params.n_quantest2
          
    output:
      set val(id), \
          val(fam), \
          val(n_quantest2), \
          val(regtrim_tree), \
          file("${id}.informative${n_quantest2}.with.${regtrim_tree}.tree*") \
          into predictions

    shell:
      """
      python3 $baseDir/bin/get_refANDinformative_seqs.py \
      --msa ${msa} --seq ${seqs} --names ${refnames} --tree ${tree} --n ${n_quantest2} \
      >> ${id}.informative${n_quantest2}.with.${regtrim_tree}.tree.aln




      quantest2 ${id}.informative${n_quantest2}.with.${regtrim_tree}.tree.aln ${ss} 

      awk 'NR==4{print}' ${id}.informative${n_quantest2}.with.${regtrim_tree}.tree.quantest2_log | cut -f2 \
      > ${id}.informative${n_quantest2}.with.${regtrim_tree}.tree.quantest2
      
      """

}


workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}" 
}