# transform FASTA files in S2SNet input files (PDB  Chain   Seq)
# @muntisa 2016

# input type:
##> CHEMBL1907597 [P47870] Gamma-aminobutyric acid
##MWRVRKRGYFGIWSFPLIIAAVCAQSVN

def Fasta2S2SNet(inFASTA,outS2SNet):
    # read the FASTA file
    finFASTA = open(inFASTA,"r")
    linesFASTA = finFASTA.readlines()
    finFASTA.close()
    Seq=""
    foutFile = open(outS2SNet,"w") # start (over)write the output file
    iline=0
    for sline in linesFASTA:
        iline+=1
        if sline[0]=='>':
            Temp = sline[2:].split(" ")
            PDBfasta=Temp[0]
            ChainFasta=Temp[1].replace("[","").replace("]","")
            if iline!=1:
                Seq=Seq+"\n"
            Seq=Seq+str(PDBfasta)+"\t"+str(ChainFasta)+"\t"
        else:
            Seq=Seq+sline[:-1]
            #print Seq
    foutFile.write(Seq)
    foutFile.close()
    return



#############################################
# main
#############################################

inFASTA="chembl_x.fa"
outS2SNet="seqs.s2snet.txt"

Fasta2S2SNet(inFASTA,outS2SNet)

print "Done! "+inFASTA+" was transformed in "+outS2SNet+"."
