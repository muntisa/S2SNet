# Transform FASTA files in S2SNet input files (PDB  Chain   Seq)
# @muntisa 2019

#input type:
##>sp|O08537.3|ESR2_MOUSE RecName: Full=Estrogen receptor beta;
##MEIKNSPSSLTSPASYNCSQSILPLEHGPIYIPSSYVESRHEYSAMTFYSPAVMNYSVPSSTGNLEGGPV

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
            if sline[-1]=="\n" : sline = sline[:-1] # removing strange characters
            if sline[-1]=="\r" : sline = sline[:-1]
            Curr = sline.split('|')
            PDBfasta=Curr[1] #sline[1:5]
            ChainFasta=Curr[2] #sline[6]
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

inFASTA="sequence.fasta"
outS2SNet="sequence.s2snet.txt"

Fasta2S2SNet(inFASTA,outS2SNet)

print ("Done! "+inFASTA+" was transformed in "+outS2SNet+".")
