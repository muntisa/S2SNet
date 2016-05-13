# Function of S2SNet
# =======================
# Authors:
#  Cristian R Munteanu:    muntisa@gmail.com
#  Humberto González-Díaz: gonzalezdiazh@yahoo.es

# input: PDB, Chain, Seq (TAB separated)

# case sensitive for the characthers
# verify if characters in Seqs can be found in Groups
# if extra chars in sequence are founded, the application is stopped
# if there are extra characters in groups there are omitted

import time
import datetime
import os, sys

from numpy import *
import math

#
# FUNCTION definitions
#

def Seq2Dot(PDB,Seq,p,g,nse,sDot): # sequence, list with position of groups, Dot file name
    if PDB=='minNetAll' or PDB=='maxNetAll' or PDB=='avgNetAll': # if are statistics files -> no labels
        sLab=' label=\"\",'
    else:
        sLab=''
    fDot= open(sDot,"w")
    colors=['blue','brown','darkorchid4','deeppink4','forestgreen','gray41','indigo','hotpink4', \
            'lightsalmon4','orangered3','orange4','royalblue4','salmon4','violetred4','yellowgreen',\
            'mediumvioletred','lightpink4','lightgoldenrod3','indianred4','chocolate1','bisque4',\
            'goldenrod4','lightslategrey','olivedrab','palevioletred4','red4','khaki3','gray17', \
            'grey','goldenrod3','darkolivegreen','burlywood']
    sDotLines=""
    sDotLines+="digraph \"net"+PDB+"\" {"
    sDotLines+="\n\tgraph [fontname = \"Helvetica\","
    sDotLines+="\n\t\tfontsize = 14,"
    #sDotLines+="\n\t\tlabel = \"\\n\\n\\n\\n"+PDB+"\\nS2SNet v3.0 (C) 2007\","
    sDotLines+="\n\t\tsize = \"7,7\", ratio=fill, rankdir=LR, center=1 ];"
    sDotLines+="\n\tnode [ color = blue, style = filled," +sLab
    sDotLines+="\n\t\tperipheries=2, fontcolor= white, fontname = \"Helvetica\" ];"
    
    fDot.write(sDotLines)

    for ip in range(len(p)):
        t=p[ip]
        for gg in range(len(g)):
            oneG=g[gg]
            if Seq[int(t[1])] == oneG[0]: # 1 because 0 is the origin (I must change this criteria; it is working with one character group 100%)
                ip=gg
            else:
                pass
        if ip <= len(colors):
            color=colors[ip]
        else:
            color='red2'
        fDot.write("\n\tnode [ color = "+color+", style = filled,")
        fDot.write("\n\t\tperipheries=2, fontcolor= white, fontname = \"Helvetica\" ];")
        
        for iip in range(len(t)-1):
            fDot.write("\n\t\""+Seq[int(t[iip])]+str(t[iip])+"\" -> \""+Seq[int(t[iip+1])]+str(t[iip+1])+"\";")
    if nse==1: #if embeded
        for iSeq in range(len(Seq)-1):
            fDot.write("\n\t\""+Seq[iSeq]+str(iSeq)+"\" -> \""+Seq[iSeq+1]+str(iSeq+1)+"\" [ color=\"red\" ];")
    fDot.write("\n\t\"00\" [ color = black ];")
    fDot.write("\n}")       
    fDot.close()
    return

def Dot2Figure(sEXE, sType, sDot, sFig, o_dir):
    cur_dir=os.getcwd()
    os.chdir(o_dir)
    cmd = '%s.exe -T%s %s -o %s' % (sEXE, sType, sDot, sFig ) # plotting graphs
    try:
        os.system(cmd)
    except:
        dlg = wx.MessageDialog(self, 'Error creating the file: "'+sFig+'". Please try again.','S2SNet: PNG File Error', wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        
    os.chdir(cur_dir)
    return

def ParseSeq(Seq,np):
    ListSubSeq=[]
    for i in range(len(Seq)-np+1):
        temp=[]
        for j in range(np):
            temp.append(Seq[i+j])
        ListSubSeq.append(temp)
    return ListSubSeq

def unique(s):
    n = len(s)
    if n == 0:
        return []
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]
    # Brute force is all that's left.
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u

def verifySeqs(sFile, gFile): #verify the characthers in Seqs with groups
    # read Group file
    fgFile = open(gFile,"r")
    glines = fgFile.readlines()
    fgFile.close()

    k=len(glines) # no of groups in input file
    g="" # Get the groups
    for gline in glines:
        if gline[-1]=="\n" : gline = gline[:-1] # removing strange characters
        if gline[-1]=="\r" : gline = gline[:-1]
        g+=gline

    # read Seqs file 
    fsFile = open(sFile,"r")
    slines = fsFile.readlines()
    fsFile.close()

    np=1
    nSeq=0
    errors={}
    err=0
    for sline in slines : # for each Seq make all the calculations
        if sline[-1]=="\n" : sline = sline[:-1] # removing strange characters
        if sline[-1]=="\r" : sline = sline[:-1]
        try:
            (PDB, Chain, Seq) = sline.split("\t") # code, chain, seq
        except:
            print "\nThe input file should have 3 TAB-delimited columns: first label, second labels and the sequence. Please correct your data. The current sequence will not be processed."
            return 1
        nSeq+=1
        uS=unique(ParseSeq(Seq,np))
        uG=unique(ParseSeq(g,np))
        extra=""
        xx=1 # for all errors
        for i in range(len(uS)):
            x=0 #for one character error
            for j in range(len(uG)):
                if uS[i]==uG[j]:
                    x=1
                    continue
            if x==0:
                extra+=str(uS[i])
                xx=0
        if xx==0:
            errors[err]=(PDB,Chain,extra)
            err+=1
    if len(errors)!=0:
        errFile="errors.txt"
        ferrFile = open(errFile,"w")
        print "Error!"
        print "(characters in Sequence but not in Groups)"
        for xe in range(len(errors)):
            print errors[xe][0],errors[xe][1],errors[xe][2]
            ferrFile.write(errors[xe][0]+" "+errors[xe][1]+" "+errors[xe][2]+"\n") 
        print "\nPlease correct the Sequence file and try again."
        print "\nYou can find the errors in ERRORS.TXT file."
        ferrFile.close()
        #sExit=raw_input('\n ... press any key to exit ...')
        return 1
    return 0

def printM(M):
    strM=''
    l=M.shape[0]
    c=M.shape[1]
    for ll in range(l):
        for cc in range(c):
            if cc==c-1:
                strM+=" "+str(M[ll][cc])+"\n"
            else:
                strM+=" "+str(M[ll][cc])
    return strM

def printV(M):
    strM=''
    l=len(M)
    for ll in range(l):
        if ll==ll-1:
            strM+=" "+str(M[ll])+"\n"
        else:
            strM+=" "+str(M[ll])
    return strM

def S2SG(Seq,k,g): # sequence to graph, k=no of groups, g=list of groups
    n=len(Seq)     # sequence length
    ip=0
    p={} # position vectors list initialization
    statGr=zeros((k),dtype=float32)
    for kk in range(k): # for each group it is a positions copy
        pp=0
        temp={}         # temp for a position group
        temp[0]=0
        l=len(g[kk])    # each group length
        for c in range(n):             # for each character in Seq
            for gg in range(l):        # for each item in a group
                if Seq[c]==g[kk][gg]:  # letter comparison of Seq with a group
                    pp+=1
                    temp[pp]=c+1       # keep the index in Seq
        statGr[kk]=len(temp)-1 #excluding the origin because is repeted in all temp[k]
        if len(temp)>1: #if the chain of the star graphs contains more than 0
            p[ip]=temp
            ip+=1
    n=n+1                     # length Seq + ZERO
    M = zeros((n,n),dtype=float32)   # initialize M with 0s
    k=len(p)                  # p position list length
    for kk in range(k):       # for each p
        temp={}
        temp=p[kk]            # working with temp=p
        pp=len(temp)          # length of current temp(p)
        for npp in range(pp):
            if npp < pp-1 :
                nl=temp[npp]
                nc=temp[npp+1]
                M[nl][nc]=1
                M[nc][nl]=1
    return (M,p,statGr)

def MatEm(M): # add embeded info to original M
    nM=M.shape[0]
    for s in range(1,nM):
        if s < nM-1:
            l=s
            c=s+1
            M[l][c]=1
            M[c][l]=1
    return M

def SWeight_old(Seq,M,w): # weights only in diagonal; Mw=weight contribution and modified M
    Seq='-'+Seq # add "-" at the begining of the Seq in order to have the same dim with M
    nM=M.shape[0]
    Mw=zeros((nM,nM),dtype=float32) # only weight in diagonal
    for iM in range(1,nM): # completeaza diagonala cu greutatea puncte
        for iw in range(1,len(w[0])):
            if Seq[iM]==w[0][iw]:
                M[iM][iM]=w[1][iw]
                Mw[iM][iM]=w[1][iw]
    return (M,Mw)

def SWeight(Seq,M,w): # Mw=weight contribution and modified M
    Seq='-'+Seq # add "-" at the begining of the Seq in order to have the same dim with M
    nM=M.shape[0]
    Mw=zeros((nM,nM),dtype=float32) # only weight in diagonal

    for iM in range(0,nM): # create Mw with diagonal as weights
        for iw in range(1,len(w[0])):
            if Seq[iM]==w[0][iw]:
                Mw[iM][iM]=w[1][iw]
                M[iM][iM]=w[1][iw]

    for i in range(nM): # modificar M with weigths M[i][j] replaced with Mw[j][j] (not only diagonal as before)
        for j in range(nM):
            if M[i][j]==1.0:
                M[i][j]=1+Mw[j][j]
    return (M,Mw)

def MPower(M,np): # power matrix
    nM=M.shape[0]
    if np==0: # zero power
        Mfin=identity(nM) # identity matrix
    if np==1:
        Mfin=M.copy()
    if np>1:  # non-zero power
        Mfin=matrix(M.copy())**int(np) # matrix to power k
##        Mfin=M.copy()
##        for npi in range(1,np):
##            Mfin=matrixmultiply(Mfin,M) # matrix multiplication
    return Mfin

def DistM(Seq,p,nse,nwi): # distances matrix
    Seq='0'+Seq # includes origin in Seq
    l=len(Seq)
    d=zeros((l,l),dtype=float32) # initialization of distances matrix
    if nse==1: #embeded ON
        for i in range(l):
            for j in range(l):
                d[i][j]=abs(i-j)

    if nse==0: # embeded OFF
        dg=zeros((l,l),dtype=float32) # distances matrix
        for i in range(l):
            for j in range(l):
                d[i][j]=abs(i-j) # distances in Seq
        for kk in range(len(p)):
            for ppi in range(len(p[kk])):
                for ppj in range(len(p[kk])):
                    dg[p[kk][ppi]][p[kk][ppj]]=float(abs(ppi-ppj)) # distance matrix in groups
        d=d+dg+2 # sum distance matrices
    return d

def DegreeM(Seq,p,nse,nwi,Mconect): # degree matrix !!!! eroare in origine !!!!
    Seq='0'+Seq # includes origin in Seq
    l=len(Seq)
    deg=zeros((l),dtype=float32) # degree matrix
    for i in range(l): # deg(i) is the i row elements sum
        sumDeg=0
        for j in range(l):
            sumDeg+=Mconect[i][j]
        deg[i]=sumDeg    
    return deg

def PairDistM(d,np): # pairs distance matrix, d=distance matrix, np=power matrix
    l=d.shape[0]
    dp=zeros((l,l),dtype=float32) # distances matrix
    for i in range(l):
        for j in range(l):
            if d[i][j]==np:
                dp[i][j]=1.0
    return dp

def HararyNo (d,nwi,Mw): # 1/2 sum[(1/dij)*wj]
    nH=0.0
    l=d.shape[0]
    for i in range(l):
        for j in range(i+1,l): # it is OK because if i=j, dij=0
            if d[i][j] != 0:
                nH+=Mw[j][j]/float(d[i][j])
    return nH

def WienerIndex (d,nwi,Mw): # 1/2 sum[(dij)*wj]
    nW=0.0
    l=d.shape[0]
    for i in range(l):
        for j in range(i+1,l):
            if d[i][j] != 0: # can be removed ...
                nW+=Mw[j][j]*d[i][j]
    return nW

def GutmanIndex(d,deg,nwi,Mw): # sum(degi*degj*Mwj/dij)deg=degree, d=dist, Mw=weights    
    nG=0.0
    l=len(deg)
    for i in range(l):
        for j in range(l):
            if d[i][j] != 0:
                nG+=Mw[j][j]*float(deg[i]*deg[j])/d[i][j]
    return nG

def SchultzIndex(d,deg,nwi,Mw): # 1/2 sum( (degi+degj)*dij*Mwj ) deg=degree, d=dist    
    nS=0.0
    l=len(deg)
    for i in range(l):
        for j in range(i+1,l):
            nS+=Mw[j][j]*float(deg[i]+deg[j])*d[i][j]
    return nS

def ATS(d,np,Mw): # ATSd= sum( dij_np*Mwi*Mwj), d=dist, np=power, Mw=weights ! only with weight
    dp=PairDistM(d,np) #new distances matrix with 1/0, where the distance is np, not only 1
    nATS=0.0
    l=dp.shape[0]
    for i in range(l):
        for j in range(l):
            if dp[i][j] != 0.00:
                nATS+=dp[i][j]*Mw[i][i]*Mw[j][j]
    return nATS/float(2.0) # only connected nodes => 12 =21

def Balaban(Mconect,d,nwi,Mw): # (n+1)sum[(aij*si*sj)*Mwj] i,j bonded, si=sum(dik)=node dist degree; # d=dist, Mw=weights, ! only with weight (without cycles no)
    l=d.shape[0]
    nBal=0.0
    edges=0.0
    for i in range(l):
        for j in range(i+1,l): # it is Ok because if i=j, i and j can not be bonded
            if Mconect[i][j]==1.0: #if i and j are bonded; includes embeded
                edges+=1
                sumI=0.0
                sumJ=0.0
                for k in range(l):
                    sumI+=d[i][k]
                    sumJ+=d[k][j]
                if sumI!=0.0 and sumJ!=0.0:
                    if nwi==1:
                        nBal+=Mw[j][j]*(sumI*sumJ)**(0.5)
                    else:
                        nBal+=(sumI*sumJ)**(-0.5)
    nBal=nBal*edges/(edges-l+2) #edges-l+1=the cyclomatic number of the molecular graph
    return nBal

def KierHall(Mconect,deg,npi,nwi,Mw,p): # d=dist, Mw=weights, ! only with weight (without cycles no),npi= a value of np
    nKH=0.0
    l=Mconect.shape[0]
    if npi==0: # power=0
        for i in range(len(deg)):
            nKH+=Mw[i][i]/math.sqrt(deg[i])
    if npi==1: # Randic index, power=1
        for i in range(l):
            for j in range(i+1,l):
                if Mconect[i][j]==1.0: #if i and j are bonded; includes embeded
                    nKH+=Mw[j][j]/math.sqrt(deg[i]*deg[j])
    if npi==2:
        for i in range(l):
            for j in range(i+1,l):
                if Mconect[i][j]==1.0:
                    for k in range(j+1,l):
                        if Mconect[j][k]==1.0:
                            nKH+=Mw[k][k]/math.sqrt(deg[i]*deg[j]*deg[k])
    if npi==3:
        for i in range(l):
            for j in range(i+1,l):
                if Mconect[i][j]==1.0:
                    for k in range(j+1,l):
                        if Mconect[j][k]==1.0:
                            for m in range(k+1,l):
                                if Mconect[k][m]==1.0:
                                    nKH+=Mw[m][m]/math.sqrt(deg[i]*deg[j]*deg[k]*deg[m])
    if npi==4:
        for i in range(l):
            for j in range(i+1,l):
                if Mconect[i][j]==1.0:
                    for k in range(j+1,l):
                        if Mconect[j][k]==1.0:
                            for m in range(k+1,l):
                                if Mconect[k][m]==1.0:
                                    for o in range(m+1,l):
                                        if Mconect[m][o]==1.0:
                                            nKH+=Mw[o][o]/math.sqrt(deg[i]*deg[j]*deg[k]*deg[m]*deg[o])                
    if npi==5:
        for i in range(l):
            for j in range(i+1,l):
                if Mconect[i][j]==1.0:
                    for k in range(j+1,l):
                        if Mconect[j][k]==1.0:
                            for m in range(k+1,l):
                                if Mconect[k][m]==1.0:
                                    for o in range(m+1,l):
                                        if Mconect[m][o]==1.0:
                                            for q in range(o+1,l):
                                                if Mconect[o][q]==1.0:
                                                    nKH+=Mw[q][q]/math.sqrt(deg[i]*deg[j]*deg[k]*deg[m]*deg[o]*deg[q])
    return nKH

def MatrixSh_old(M,deg,Mw): # Shannon Entropy Index ... -sum(pi*log(pi))
    Sh=0.0
    n=M.shape[0]
    V=ones((n),dtype=float32)/float(n)
    P=matrixmultiply(V,M) # dot product of matrix M and vector V (matrixmultiply is faster than DOT)
    #in the previouse version the product was M,V
    
    sumSh=0
    for i in range(len(P)):
        sumSh+=P[i]*log(P[i])
    Sh=-sumSh

    return Sh

def MatrixSh(M,nwi,Mw,deg): # Shannon Entropy Index ... -sum(pi*log(pi))
    # Sh[npi]=MatrixSh(MpN[npi],nwi,Mw,deg)
    Sh=float(0)
    n=M.shape[0]
    #unique V = 1/n vector
    V=ones((n),dtype=float32) ## V=ones([n], Float)
    P=zeros((n),dtype=float32)
    #in the case of no weigths Mw is full of ones
    V=diagonal(Mw)/float(trace(Mw))
    # eliminated in order to have the same 1/n vector V
    
##    if nwi==1:
##        V=diagonal(Mw)/float(trace(Mw))
##    else:
##        #V=diagonal(Mw)/float(trace(Mw))
##        sumDeg=0
##        for i in range(len(deg)):
##            sumDeg+=deg[i]
##        for i in range(len(deg)):
##            if sumDeg!=0:
##                V[i]=float(deg[i])/float(sumDeg)
##            else:
##                V[i]=0
    
    P=dot(V.copy(),M.copy()) # P=matrixmultiply(M,V) # dot product of matrix M and vector V (matrixmultiply is faster than DOT)
    if P.shape[0]==1:
        P.resize(n)
    # in the next version we need to chack the order of the multiplication!!!!!!!!!!! old version = M,V; new vers=V,M
    sumSh=0
    for i in range(len(P)):
        if P[i]!=0:
            sumSh+=P[i]*log(P[i])
    Sh=-sumSh
    return Sh

def S2SGmain(emb, wts, dot, markov, det, power, sequences, groups, weights,results,details,orig_dir):
    # MAIN
    tt1=time.clock()

    #print the parameters of the calculations ...........
    
    # INPUT
    sFile = sequences
    gFile = groups
    
    print "\n* Verifying the characters in Sequences and Groups ... ",
    iVer=verifySeqs(sFile,gFile) #compare Seqs with groups
    # if Verification is 1 = exist errors! 
    if iVer==1:
        return # exit the calculation
    print "OK!"

    # read Seqs file
    fsFile = open(sFile,"r")
    slines = fsFile.readlines()
    fsFile.close()
    
    # read Group file
    fgFile = open(gFile,"r")
    glines = fgFile.readlines()
    fgFile.close()

    k=len(glines) # no of groups in input file

    g={} # Get the groups
    gi=0
    for gline in glines:
        if gline[-1]=="\n" : gline = gline[:-1] # removing strange characters
        if gline[-1]=="\r" : gline = gline[:-1]
        g[gi]=gline
        gi+=1

    # POWER of M
    np=int(power)
    #
    # EMBEDED
    if emb==True:
        nse=1
    else:
        nse=0
    # 
    # WEIGHTs
    #
    if wts==True:
        nwi=1
        wFile=weights
        # read Weights file
        fwFile = open(wFile,"r")
        wlines = fwFile.readlines()
        fwFile.close()
        lW=len(wlines)
        w2 = zeros((lW+1),dtype=float32) # weight values
        w1={}
        w2[0]=0.00
        w1[0]="-"
        ww=1
        for wline in wlines :
            if wline[-1]=="\n" : wline = wline[:-1] # removing strange characters
            if wline[-1]=="\r" : wline = wline[:-1]
            (sl,sw) = wline.split("\t") # letter, weigth
            w2[ww]=float(sw)
            w1[ww]=sl
            ww+=1
        w=[w1,w2]
    else:
        nwi=0
        Mw={} # Null matrix of weights even if weight is OFF

    #
    # Normalization type:
    # 0 = non-Markov = M**np -> normalization of each -> trace of each
    # 1 = possible Markov = normalization of M -> power np -> trace of each
    #
    if markov==True:
        Mark=1
    else:
        Mark=0

    #
    # RESULTS file
    #
    rFile = results
    if det == True:
        nDet=1
        dFile = details
    else:
        nDet=0

    # --------------------------------------------------------------------------------
    #  PROCESSING
    # --------------------------------------------------------------------------------

    # open Summary results file (overwrite)        
    frFile = open(rFile,"w")

    if nDet==1:
        # open Details results file (overwrite)
        fdFile = open(dFile,"w")

    # dinamic table header
    sHeader="PDB\tChain\tSeq" #modified !
    if Mark==1:
        for i in range(0,np+1):
            sHeader+="\tSh"+str(i)
    for i in range(0,np+1):
        sHeader+="\tTr"+str(i)
    sHeader+="\tH\tW\tS6\tS"
    if nwi==1:
        for i in range(1,np+1):
            sHeader+="\tATS"+str(i)
    sHeader+="\tJ"
    for i in range(np+1): # limited to 0-power
        if i==1:
            sHeader+="\tX"+str(i)+"(R)"
        else:
            sHeader+="\tX"+str(i)
    sHeader+="\n" # return new line
    frFile.write(sHeader)


    print "* Processing:\n"
    nSeq=0
    SeqList=[]
    allSeq=len(slines)
    
    # initialize entire list for min, max and avg graph sequences    
    statGrAll={}
    for sline in slines : # for each Seq make all the calculations
        if sline[-1]=="\n" : sline = sline[:-1] # removing strange characters
        if sline[-1]=="\r" : sline = sline[:-1]
        try:
            (PDB, Chain, Seq) = sline.split("\t") # code, chain, seq
            if len(Seq)>5000:
                print "The "+str(PDB)+" is too large for this application and will be skipped"
                continue
        #except ValueError:
        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise
            #print "The input file should have 3 TAB-delimited columns: first label, second labels and the sequence. Please correct your data. The current sequence will not be processed."
            #continue
        nSeq+=1
        SeqList.append(PDB+Chain)
        sOut= "\n%d of %d (%d%%): %s%s" % (nSeq, allSeq,nSeq*100/allSeq,PDB,Chain)
        if nDet==1:
            fdFile.write(sOut+"\nSeq="+Seq+"\n") 
        print sOut
        sOut= str(PDB)+"\t"+str(Chain)+"\t"+Seq
        frFile.write(sOut,)

        nM=len(Seq)+1 # add the origin

        # Seq to matrix
        (M,p,statGr)=S2SG(Seq,k,g)
        statGrAll[nSeq-1]=statGr
        
        if dot==True:
            # only DOT file
            sDot=orig_dir+"dot\\"+str(PDB)+str(Chain)+".txt"
            try:
                fDot= open(sDot,"w")
                fDot.close()
            except:
                print "The labels for "+str(PDB)+" contain strange characters and the DOT and PNG files can not be created. Please correct your input file and restart the calculation."
                continue

            Seq2Dot(str(PDB)+str(Chain),"0"+Seq,p,g,nse,sDot)
            
        if nDet==1:
            fdFile.write("\nM=conectivities matrix\n")
            fdFile.write(printM(M))

        Mconect=zeros((nM,nM),dtype=float32) # initialize Mconect
        if nse==0:
            for l in range(nM):
                for c in range(nM):
                    Mconect[l][c]=M[l][c]
        # check embeded
        if nse == 1:
            M=MatEm(M) 
            for l in range(nM):
                for c in range(nM):
                    Mconect[l][c]=M[l][c]
            if nDet==1:
                fdFile.write("\nMemb=embeded modified matrix M\n")
                fdFile.write(printM(M))
                
        Mw=zeros((nM,nM),dtype=float32)
        # check weight
        if nwi == 1:
            M=SWeight(Seq,M,w)[0]
            Mw=SWeight(Seq,M,w)[1]
            if nDet==1:
                fdFile.write("\nMweights=weights modified matrix M\n")
                fdFile.write(printM(M))
                fdFile.write("\nMw=weights only matrix M\n")
                fdFile.write(printM(Mw))
                
        # Normalization type
        if Mark==0:
            # generate power of M
            Mp={} # list of matrices M^i
            for npi in range(0,np+1):
                Mp[npi]=MPower(M,npi)
                if nDet==1:
                    fdFile.write("\nMpower["+str(npi)+"]\n")
                    fdFile.write(printM(Mp[npi]))
            # matrix normalize
            if nDet==1:
                fdFile.write("\n* Normalize matrices\n")
            nMp=len(Mp)
            MpN={} # list of normalized matrices M^i
            for nMpi in range(nMp):
                MM=(Mp[nMpi]).copy()
                for li in range(nM):
                    tli=0.0
                    for co in range(nM):
                        tli+=MM[li][co]
                    for co in range(nM):
                        if tli!=0:
                            MM[li][co]=float(MM[li][co])/float(tli)
                        else:
                            MM[li][co]=0.0 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                MpN[nMpi]=MM.copy()
                if nDet==1:
                    fdFile.write("\nM["+str(nMpi)+"] normalized\n")
                    fdFile.write(printM(MpN[nMpi]))
                    
        if Mark==1:
            # matrix normalize
            MM=zeros((nM,nM),dtype=float32)
            for li in range(nM):
                tli=0.0
                for co in range(nM):
                    tli+=M[li][co]
                for co in range(nM):
                    if tli!=0:
                        MM[li][co]=float(M[li][co])/float(tli)
                    else:
                        MM[li][co]=0.0            
            if nDet==1:
                fdFile.write("\n* Normalize matrix M\n")
                fdFile.write(printM(MM))
            # generate power of M
            MpN={} # list of matrices M^i
            for npi in range(0,np+1):
                MpN[npi]=MPower(MM,npi)
                #print "trace=",trace(MpN[npi])
                if nDet==1:
                    fdFile.write("\nMpower["+str(npi)+"]\n")
                    fdFile.write(printM(MpN[npi]))
        
        ###########################################################
        #
        # QSAR numbers
        #
        ###########################################################
     
        if nDet==1:
            fdFile.write("\n* Topological Indices\n")
            
        d=DistM(Seq,p,nse,nwi)               # distance matrix
        deg=DegreeM(Seq,p,nse,nwi,Mconect)   # degree matrix

        if nDet==1:
            fdFile.write("\n\nd = distance matrix\n")
            fdFile.write(printM(d))
            fdFile.write("\ndeg = degree matrix\n")
            fdFile.write(printV(deg))
            fdFile.write("\n\nMconect=connectivity matrix depending of embeded\n")
            fdFile.write(printM(Mconect))
            
        if nwi==0:
            Mw=ones((nM,nM),dtype=float32)
            
        # if Markov, calculate Shannon Entropy indices
        if Mark==1:
            Sh={}
            Sh[0]=0.0
            for npi in range(0,np+1):
                Sh[npi]=MatrixSh(MpN[npi],nwi,Mw,deg)
                frFile.write("\t"+str(Sh[npi]),)
                if nDet==1:
                    fdFile.write("\nSh["+str(npi)+"]="+str(Sh[npi]))
            if nDet==1:
                fdFile.write("\n")
                
        # trace (nTpow)
        nT={}
        for iMpN in range(len(MpN)):
            nT[iMpN]=trace(MpN[iMpN])
            if nDet==1:
                fdFile.write("\nTr["+str(iMpN)+"]="+str(nT[iMpN]))
            frFile.write("\t"+str(nT[iMpN]),)
                  
        # Harary number (nH)
        nH=HararyNo(d,nwi,Mw)
        if nDet==1:
            fdFile.write("\n\nH="+str(nH))
        frFile.write("\t"+str(nH),)

        # Wiener index (nW)
        nW=WienerIndex(d,nwi,Mw)
        if nDet==1:
            fdFile.write("\n\nW="+str(nW))
        frFile.write("\t"+str(nW),)

        # Gutman Index (nG)
        nG=GutmanIndex(d,deg,nwi,Mw)
        if nDet==1:
            fdFile.write("\n\nG="+str(nG))
        frFile.write("\t"+str(nG),)

        # Schultz Index (nS)
        nS=SchultzIndex(d,deg,nwi,Mw)
        if nDet==1:
            fdFile.write("\n\nS="+str(nS)+"\n")
        frFile.write("\t"+str(nS),)

        # Moreau-Broto ATS
        if nwi==1: # only with Weight ON 
            for iATS in range(1,np+1):
                dp=PairDistM(d,iATS)# pair matrix
                nATS=ATS(d,iATS,Mw)
                if nDet==1:
                    fdFile.write("\nATS"+str(iATS)+"="+str(nATS))
                frFile.write("\t"+str(nATS),)

        # Balaban
        if nDet==1:
            fdFile.write("\n\n")
        nBal=Balaban(Mconect,d,nwi,Mw)
        if nDet==1:
            fdFile.write("J="+str(nBal))
        frFile.write("\t"+str(nBal),)

        # Kier-Hall
        if nDet==1:
            fdFile.write("\n")
        for npi in range(np+1): # only 0,1,2 indexes (1+1)
            nKH=KierHall(Mconect,deg,npi,nwi,Mw,p)
            if npi==1: #for Randic index only
                if nDet==1:
                    fdFile.write("\nX"+str(npi)+" (Randic) = "+str(nKH))
            else:
                if nDet==1:
                    fdFile.write("\nX"+str(npi)+"="+str(nKH))
            frFile.write("\t"+str(nKH),)
        frFile.write("\n") # end for one Seq
        if nDet==1:
            fdFile.write("\n---------------------------------------------------------------")
    frFile.close() # close results file
    if nDet==1:
        fdFile.close()

    # list the total execution time
    tt2=time.clock()
    print "\nTotal Execution Time: %(ddiffsec).1fs\n\n\n" %\
          {"ddiffsec": tt2-tt1}
    
    return (SeqList,statGrAll,g)
