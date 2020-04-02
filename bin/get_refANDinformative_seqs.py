#!/usr/bin/env python3

__description__ = '''

Get reference sequences + (N-ref) informative sequences for --> QuanTest2

Warning: You should have T-coffee installed before running this script.

'''

def get_informative_fasta(fasta,tree,n):

    ''' Get informative n sequences from a fasta file <msa or seq>, based on the guide tree <tree> 
    
    An ordered dictionary {name : sequence} with the n informative sequences is returned.
    '''

    import collections

    inf_fasta=collections.OrderedDict()
    
    tmp=os.popen("t_coffee -other_pg seq_reformat -in {} -in2 {} -action +regtrim {}".format(fasta,tree,n)).read()
    tmp=[i for i in tmp.split(">") if i]

    for item in tmp:
        item=[i for i in item.split("\n") if i]
        name=item[0]
        seq=""
        for i in range(1,len(item)):
            seq+=item[i]
        inf_fasta[name]=seq

    return(inf_fasta)


def merge_ref_informative_seqs(ref_fasta,inf_fasta,n):
    ''' Merge references + (n-ref) informative sequences 

    An ordered dictionary {name : sequence} with n sequences (ref + informative) is returned.
    '''

    import collections

    # Reference sequences
    out=ref_fasta.copy()
    # Add N-ref informative sequences
    n=len(ref_fasta)
    for name,seq in inf_fasta.items():
        if n>=args.n:
            break
        if name in ref_fasta.keys():
            continue
        out[name]=seq
        n+=1
    
    return(out)



if __name__ == '__main__':

    import argparse
    import os
    import sys
    from get_seqs import get_seqs
    from readFilesFunctions import read_fasta, read_aux
   

    app = argparse.ArgumentParser(description=__description__)
    app.add_argument('--fasta',type=str,help="Fasta file <seq or msa>")
    app.add_argument('--names',type=str,help="A file with a list of reference sequence names. \
                     In particular, a total of 3 sequences with known structures are assigned as the references \
                     (whose secondary structure will be predicted using QuanTest2)")
    app.add_argument('--tree',type=str,help="Guide tree to be used to retrieve the informative sequences --> T-coffee +regtrim")
    app.add_argument('--n',type=int,help="Total number of sequences (ref + informative)")
    args = app.parse_args() 

    # Read files
    fasta=read_fasta(args.fasta)
    names=read_aux(args.names)
    
    # Get reference sequences
    ref_fasta=get_seqs(fasta,names)
    # Get N-ref informative sequences
    inf_fasta=get_informative_fasta(args.fasta,args.tree,args.n)
    # Merge ref + informative sequences
    out=merge_ref_informative_seqs(ref_fasta,inf_fasta,args.n)

    # Write output: references + (N-ref) informative sequences
    for name,seq in out.items():
        sys.stdout.write(">"+name+"\n"+seq+"\n")

