#!/usr/bin/env python3
__author__ = "fabian"
__date__ = "2019"


# Note: this script has been modified in order to use a standalone version of PSIPRED to predict secondary structures (SJ March 2020)
# Usage: python3 quantest2Local.py <alignment file> <refs ss file>


# Accessory scripts from other programs
# Please set the path before running the script
reformatPl="/usr/src/app/hhsuite2/scripts/reformat.pl"      # HH-suite reformat.pl script
a3m2mtxPl="/usr/src/app/hhsuite2/scripts/a3m2mtx.pl"        # DeepMSA: a3m2mtx.pl script
runpsipredSh="/usr/src/app/psipred/runpsipred"              # Script to run PSIPRED


GAPCHARS = ['-', '.']
GAPTHRESH = 0.90
REF_TARGET = 3 # exepect every reference file to have that many references
verbose =  True  # minimum verbosity
VERBOSE =  False #True ## default should be 'False'
VVERBOSE = False #True ## default *definitely* should be 'False'
RUTHLESS = True #False #True ## tidy up files down-loaded from Jpred, should be 'True'
UNKNOWN = -1
import sys
import subprocess

nFiles_2 = 0
alnFiles = []
ssFiles  = []
cnt = 0
indiScores = []
predicted = []

myBasename=sys.argv[0].rsplit('/',1)[0]
if VERBOSE:
    print("myBasename = {}".format(myBasename))



    

#### Commandline() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# put commandline arguments into proper arrays/variables
# return nFiles_2, number of alignment files (= number of secondary structure files)
def Commandline(argv):

    usage = "Usage: "+argv[0]+" <alignment-file1> ... <secondary-file1> ... "
    
    argaux = []
    
    for a in argv[1:]:
            argaux.append(a)  # Arguments: <alignment-file1> ... <secondary-file1> ...

    if VVERBOSE:
        print(argaux)
    
    nFiles = len(argaux)
    nFiles_2 = int(nFiles / 2)
    if nFiles < 2:
        print("Commadline Error: at least 2 input files expected")
        print(usage)
        quit()
    if nFiles % 2 != 0:
        print("equal number of alignment files and secondary structure files expected")
        print(usage)
        quit()
    # 1st half of input files are alignments, 2nd half of input files should be secondary structures
    for i in range(nFiles_2):
        alnFiles.append(argaux[i])
        ssFiles.append(argaux[i+nFiles_2])

    return nFiles_2
# end of def Comanndline() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#### ReadSS() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# read file with secondary structure information
# return NOOFREF, number of secondary structures
def ReadSS(ssFilePtr, refLabel, refSeq):
    
    NOOFREF = UNKNOWN
    for line_ in ssFilePtr:
        line = line_.rstrip()
        if line == "":
            continue
        if line[0] == ">":
            NOOFREF = NOOFREF+1
            refLabel.append(line[1:])
            refSeq.append("")
        else:
            refSeq[NOOFREF] =  refSeq[NOOFREF]+line
    ssFilePtr.close()
    NOOFREF = NOOFREF+1

    if VVERBOSE:
        print(refSeq)

    return NOOFREF


# end of ReadSS() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#### def ReadAln() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# read alignment file
# return number of sequences in alignment file
def ReadAln(alnFilePtr, label, seq, cnt, NOOFREF, refLabel, REFPOS, family):

    cnt = UNKNOWN
    for line_ in alnFilePtr:
        line = line_.rstrip()
        if line == "":
            continue
        if line[0] == ">":
            cnt = cnt+1
            label.append(line[1:])
            seq.append("")        
            for i in range(NOOFREF):
                if label[cnt] == refLabel[i]:
                    REFPOS[i] = cnt
        else:
            if cnt == UNKNOWN:
                print("Format Error: no valid sequence label found")
                quit()
            seq[cnt] = seq[cnt]+line

    cnt = cnt+1
    
    if VVERBOSE:
        for i in range(cnt):
            print("seq[{}] = {}\t{}".format(i, label[i],seq[i]))
            
    
    # check that all reference sequences were present
    for i in range(NOOFREF):
        if REFPOS[i] == UNKNOWN:
            print("WARNING: reference sequence {} ({}) not in alignment".format(i, refLabel[i]))
    
    # check all sequence lines have the same length
    alnLength = len(seq[0])
    if VVERBOSE:
        print("1st sequence is {} long, there are {} sequences".format(alnLength,cnt))
    for i in range(cnt):
        if len(seq[i]) != alnLength:
            print("WARNING: sequence {} ({}) has length {}, not same as previous sequence/s {}".format(i,label[i],len(seq[i]),alnLength))
            quit()
    if VVERBOSE:
        logFile.write("all aligned sequences in {} have correct lengths\n".format(family))
    
    return cnt

# end of ReadAln() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



#### def Gblocks() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# put sequence to be predicted on top, remove columns in alignment where 1st sequences has gaps
# write this construct to file, return name of this file
# So, given 3 reference sequences per family, 3 .blk files containing the modified alignments will be created
def Gblocks(family, ri, refLabel, refSeq, REFPOS, cnt, seq):

    alnLength = len(seq[0])
    
    blockfile=family+"-"+str(ri)+".blk"
    try:
        blk = open(blockfile, 'w')
    except:
        print("Could not open file {} for writing (permission?)".format(blockfile))
        quit()
            
    # print sequence that should be predicted
    if VVERBOSE:
        print("ref {}: {}, position {}, write to {}".format(ri,refLabel[ri],REFPOS[ri],blockfile))
    printline = ""
    for ci in range(alnLength):
        c = seq[REFPOS[ri]][ci]
        if c != '-' and c != '.':
            printline += c
    # check that the number of residues equals number of ss states
    if len(printline) != len(refSeq[ri]):
        print("number of residues ({}) for reference sequence {} ({}[{}]) does not equal #SS states ({})".format(len(printline),refLabel[ri],family,ri,len(refSeq[ri])))
        quit()
    elif VVERBOSE:
        logFile.write("reference sequence {}/{} has correct number of residues/states\n".format(refLabel[ri],family))
    blk.write(">{}\n{}\n".format(refLabel[ri],printline))
    if VVERBOSE:
        print(printline)
    # do the other sequences
    k=1
    for j in range(cnt):
        if j == REFPOS[ri]:
            continue
        printline = ""
        for ci in range(alnLength):
            cr = seq[REFPOS[ri]][ci]
            if cr != '-' and cr != '.':
                printline += seq[j][ci]
        gapCnt = 0
        for gc in GAPCHARS:
            gapCnt = gapCnt + printline.count(gc)
        if (gapCnt / len(refSeq[ri])) > GAPTHRESH:
            continue
        blk.write(">{}\n{}\n".format(k,printline))
        k = k+1
        if VVERBOSE:
            print("{} ({})".format(printline,j))
    blk.close()

    return blockfile

# end of Gblocks() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# extract and score predictions
#### def ExtactAndScore() >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# Extract MTX file from MSA, predict secondary structure and score the accuracy
def ExtractAndScore(blkReg, nmeReg, refReg, indiScores, predicted, logFile):
    
    import os
    import subprocess
    
    for j in range(len(blkReg)):

        # Files
        blkfile=blkReg[j]                       # Modified MSA in fasta format
        pref=os.path.splitext(blkfile)[0]       # Prefix
        a3mfile=pref+".a3m"                     # Modified MSA in a3m format
        mtxfile=pref+".mtx"                     # MTX file

        # Extract MTX file from modified MSA
        subprocess.call(["perl",reformatPl,"fas","a3m",blkfile,a3mfile],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        subprocess.call(["perl",a3m2mtxPl,a3mfile,mtxfile],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

        # Predict secondary structure from MTX file
        subprocess.call(["tcsh",runpsipredSh,mtxfile],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)   # This creates .ss .ss2 and .horiz files
        ssfile=pref+".ss"

        # Read predicted ss file
        states=""
        with open(ssfile) as f:
            for line in f:
                line=line.strip("\n")
                states+=line[7]
        predicted[j]=states
        
        # Score
        match=0.0
        for c in range(len(refReg[j])):

            r = refReg[j][c]
            s = states[c]
            if   (r == 'H' or r == 'G' or r == 'I') and (s == 'H' or s == 'G' or s == 'I'):
                match = match + 1
            elif (r == 'E' or r == 'B')             and (s == 'E' or s == 'B'):
                match = match + 1
            elif (r == '-' or r == 'C')             and (s == '-' or s == 'C'):
                match = match + 1
            elif (r != 'H' and r != 'G' and r != 'I' and r != 'E' and r != 'B' and r != '-' and r != 'C') or (s != 'H' and s != 'G' and s != 'I' and s != 'E' and s != 'B' and s != '-' and s != 'C'):
                print("WARNING: unknown SS status, c={}: r = {} / s = {}".format(c,r,s))

        perc = match/len(refReg[j])*100
        indiScores[j] = perc
        if VERBOSE:
            print("there are {} matches out of {} = {}%".format(match,len(refReg[j]),perc))

        logFile.write("{}\t{}\t{}\n".format(nmeReg[j], perc, states))
    
# end of ExtractAndScore() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



def main():
    
    nFiles_2 = Commandline(sys.argv)
    if VERBOSE:
        print(alnFiles)
        print(ssFiles)
    if verbose:
        print("parsed commandline")
    
    # go through different alignment/secondary files
    eCnt=0
    blkReg = []
    nmeReg = []
    refReg = []
    
    for II in range(nFiles_2):
    
        alignment = alnFiles[II]
        aux = alignment[0:alignment.rfind(".")]
        family = aux[aux.rfind("/")+1:]
        ylimaf = family.replace('-','_').replace('.','_')

        LOGFILE = family + ".quantest2_log"
        try:
            logFile = open(LOGFILE,"w")
        except:
            print("OutputError: Could not open file {} for writing".format(LOGFILE))

        try:
            alnFilePtr = open(alignment,"r")
        except FileNotFoundError:
            print("Input Error: alignment file {} (index {}) does not exist".format(alignment,II))
            quit()
    
        secondary = ssFiles[II]
        try:
            ssFilePtr = open(secondary,"r")
        except FileNotFoundError:
            print("Input Error: alignment file {} (index {}) does not exist".format(secondary,II))
            quit()
        
        # refLabel/refSeq are different for every file, initialise inside the main loop
        refLabel = []
        refSeq = []
        NOOFREF = ReadSS(ssFilePtr, refLabel, refSeq)
        if NOOFREF != REF_TARGET:
            print("WARNING: SS-file {} has {} references, expecting {} -- possible problem".format(II, NOOFREF, REF_TARGET))

        refReg.extend(refSeq) # refReg is the list of all reference structures that is kept till the end
        
        REFPOS = []
        for i in range(NOOFREF):
            REFPOS.append(UNKNOWN)
    
    
    
        label = []
        seq = []
        alnLen = UNKNOWN
        # read fasta files, store as vienna
        cnt = UNKNOWN
        cnt = ReadAln(alnFilePtr, label, seq, cnt, NOOFREF, refLabel, REFPOS, family)
        
        if VVERBOSE:
            print("refLabel = {}".format(refLabel))  # References id
            print("refSeq = {}".format(refSeq))      # References SS
            print("#r = {}".format(NOOFREF))         # Number of references
            print("refpos = {}".format(REFPOS))      # References position in the alignment
            print("r1 = {}".format(seq[REFPOS[0]]))
            print("r2 = {}".format(seq[REFPOS[1]]))
            print("r3 = {}".format(seq[REFPOS[2]]))

        for ri in range(NOOFREF):
            # do gblocks
            blockfile = Gblocks(family, ri, refLabel, refSeq, REFPOS, cnt, seq)
            myname=ylimaf+"__"+str(ri)
            blkReg.append(blockfile)
            nmeReg.append(myname)
        
        if VVERBOSE:
            print("blk = {}".format(blkReg))
            print("nme = {}".format(nmeReg))
    
    if verbose:
        print("submitted gblocks")

    indiScores = [0]*len(blkReg)
    predicted  = [""]*len(blkReg)
    ExtractAndScore(blkReg, nmeReg, refReg, indiScores, predicted, logFile)

    if verbose:
        print("retrieved predictions")

    for j in range(len(blkReg)):
        if j % NOOFREF == 0:
            mysum = 0.00
        mysum = mysum + indiScores[j]
        if j % NOOFREF == NOOFREF-1:
            print("{}\t{}".format(alnFiles[int(j/NOOFREF)],mysum/NOOFREF))
            logFile.write("{}\t{}\n".format(alnFiles[int(j/NOOFREF)],mysum/NOOFREF))
    
    if verbose:
        print("extracted scores")


if __name__ == "__main__":
    main()
