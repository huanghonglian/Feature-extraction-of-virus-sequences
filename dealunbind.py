import pandas as pd
from glob import glob
import re
import os
import sys
import codecs
from collections import defaultdict
import re
ref= pd.read_csv('./public/reference/ref_seq.fasta', header=0,index_col=None)
ref=ref.iloc[:, 0].tolist()[0]
sample= pd.read_csv('./public/sample/metadata.txt', header=0,index_col=None)
file_name=sample["sample"].tolist()
count={} 
for i in range(1,1+len(ref)):
    count[str(i)]=0
datalast=pd.DataFrame(list(count.items()),columns=['id', 'sequence_nub'])
datalast["sample"]="sample"

if not os.path.exists('./result/unbindcleansample/'):
    os.makedirs('./result/unbindcleansample/')
if not os.path.exists('./result/unbindcleansample/datadeal/'):
    os.makedirs('./result/unbindcleansample/datadeal/')
if not os.path.exists('./result/unbindcleansample/del/'):
    os.makedirs('./result/unbindcleansample/del/')

for file in file_name:
    with open("./result/unbindsam/"+file+"_unbind.sam","r") as pdinf:
        count={} 
        for i in range(1,1+len(ref)):
            count[str(i)]=0
        repeat={}
        print ('Start processing file:',file)
        with open("./result/unbindcleansample/datadeal/"+file+".txt","w") as out:
            with open("./result/unbindcleansample/del/"+file+"int.txt","w") as out1:
                for lines in pdinf.readlines()[2:]:
                    line=lines.split("\t")
                    if line[0] not in repeat.keys():
                        repeat={}
                        repeat[line[0]]=[]
                    letraw=re.findall(r"[A-Z]",line[5])
                    start=0
                    ####S,D,DEL,I,D
                    pattern=re.compile(r"\d+[A-Z]")
                    inf=pattern.findall(line[5])
                    inset=0
                    if line[5]!="*":
                        nublast=inf[len(inf)-1]
                        seqraw=line[9]
                        if nublast[len(nublast)-1:]=="S":
                            b=len(inf)-1
                            kk=re.findall(r"\d+",inf[b])[0]
                            seqraw=seqraw[:-int(kk)]
                            inf=inf[:b]                        
                        del1=1
                        for v in range(len(inf)):
                          nub=int(re.compile(r"\d+").findall(inf[v])[0])
                          let=re.compile(r"[A-Z]").findall(inf[v])[0]
                          if let=="M":
                            inset=inset+nub
                          if let=="S" and v==0:
                            #inset=inset+nub
                            seqraw=seqraw[nub:]
                          if let=="H":
                            seqraw=seqraw
                          if let=="D":
                            inset=inset+nub
                            startseq=seqraw[0:start]
                            endseq=seqraw[start:]
                            seqraw=startseq+"-"*nub+endseq
                          if let!="I" and let!="S" and let!="H":
                            start=start+nub
                          if let=="I":
                            startseq=seqraw[:inset]
                            intseq=seqraw[inset:inset+nub]
                            insetsit=inset+int(line[3])
                            endseq=seqraw[(inset+nub):]
                            seqraw=startseq+endseq
                            out1.write(">"+line[0]+"\t"+line[5]+"\t"+str(insetsit)+"\t"+intseq+"\n")
                          if let!="D":
                            del1=del1+nub 
                    ####out result 
                    if "M" in letraw:          
                        if len(letraw)==1 and let[0]=="M":
                            dd=list(range(int(line[3]),int(line[3])+len(line[9])))
                            inter=list(set(repeat[line[0]]).intersection(set(dd)))
                            repeat[line[0]].extend(dd)
                            if len(inter)==0:
                                check="first"
                            else:
                                check="second"
                            out.write(">"+line[0]+"\t"+line[5]+"\t"+line[3]+"\t"+line[9]+"\t"+check+"\n")
                            for x in range(len(line[9])):
                                nubseq=str(int(line[3])+x)
                                count[nubseq]+=1
                        elif len(letraw)!=1:
                            if "D" in letraw or "I" in letraw or "N" in letraw or "S" in letraw or "H" in letraw: 
                                dd=list(range(int(line[3]),int(line[3])+len(seqraw)))
                                inter=list(set(repeat[line[0]]).intersection(set(dd)))
                                repeat[line[0]].extend(dd)
                                if len(inter)==0:
                                    check="first"
                                else:
                                    check="second"
                                out.write(">"+line[0]+"\t"+line[5]+"\t"+line[3]+"\t"+seqraw+"\t"+check+"\n")  
                                for x in range(len(seqraw)):
                                    nubseq=str(int(line[3])+x)
                                    count[nubseq]+=1
                # else:
                #     out1.write(lines)
        data=pd.DataFrame(list(count.items()),columns=['id', 'sequence_nub'])
        data["sample"]=file
        datalast = pd.concat((datalast,data),axis=0)
datalast=datalast.loc[datalast['sample'] != "sample"] 
datalast.to_csv('./result/datalastunbind.txt',sep='\t',index=False)
    
                       



