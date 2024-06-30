import pandas as pd
from glob import glob
import re
import os
import sys
import codecs
from collections import defaultdict
import re
sample= pd.read_csv('./public/sample/metadata.txt', header=0,index_col=None)
sample=sample["sample"].tolist()
ref= pd.read_csv('./public/reference/ref_seq.fasta', header=0,index_col=None)
ref=ref.iloc[:, 0].tolist()[0]
count={} 
for i in range(1,1+len(ref)):
    count[str(i)]=0
datalast=pd.DataFrame(list(count.items()),columns=['id', 'sequence_nub'])
datalast["sample"]="sample"

if not os.path.exists('./result/finaldata/'):
    os.makedirs('./result/finaldata/')
if not os.path.exists('./result/finaldata/seq/'):
    os.makedirs('./result/finaldata/seq/')
if not os.path.exists('./result/finaldata/count/'):
    os.makedirs('./result/finaldata/count/')
for file in sample:
    print ('Start processing file:',file)
    count={} 
    seq={}
    for i in range(1,1+len(ref)):
        count[i]=0
        seq[i]="*"
    with open("./result/finaldata/seq/"+file+".txt","w") as outfile:
        with open("./result/cleansample/"+file+".txt","r") as pdinf:
            for lines in pdinf.readlines():
                line=lines.strip("\n").split("\t")
                for x in range(len(line[3])):
                    nub=int(line[2])+x
                    count[nub]+=1
                outfile.write(line[0]+"\t"+line[2]+"\t"+line[3]+"\n")
        with open("./result/unbindcleansample/lastdata/"+file+"unbindlast.txt","r") as out:
            for lines in out.readlines():
                line=lines.strip("\n").split("\t")
                for x in range(len(line[3])):
                    nub=int(line[2])+x
                    count[nub]+=1
                outfile.write(line[0]+"\t"+line[2]+"\t"+line[3]+"\n")
        with open("./result/unbindcleansample/lastdata/"+file+"deal.txt","r") as out2:
            for lines in out2.readlines():
                line=lines.strip("\n").split("\t")
                for x in range(len(line[3])):
                    nub=int(line[2])+x
                    count[nub]+=1
                outfile.write(line[0]+"\t"+line[2]+"\t"+line[3]+"\n")
        data=pd.DataFrame(list(count.items()),columns=['id', 'sequence_nub'])
        data["sample"]=file
        datalast = pd.concat((datalast,data),axis=0)
    datalast=datalast.loc[datalast['sample'] != "sample"] 
    datalast.to_csv("./result/finaldata/count/"+file+'count.txt',sep='\t',index=False)
    
                       



