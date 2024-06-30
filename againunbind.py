
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

for file in sample:
    data= pd.read_csv("./result/unbindcleansample/datadeal/"+file+".txt", names= ["ID","match","start","seq","check"],header=None,sep="\t",index_col=None)  
    data['Total'] = data.groupby("ID").ID.transform('count')
    data['checkTotal'] = data.groupby(["ID","check"]).ID.transform('count')
    #print (data.head())
    dataunique=data.loc[data['Total'] == 1]
    dataunbind=data.loc[(data['Total'] == 2) & (data['checkTotal'] == 2)]
    dataunique = pd.concat((dataunique,dataunbind),axis=0)
    dataunique.to_csv('./result/unbindcleansample/lastdata/'+file+"unbindlast.txt",sep='\t',index=False,header=0)
    dataunique=dataunique.iloc[0:1,:]
    datadeal=data.loc[(data['Total'] == 2) & (data['checkTotal'] != 2)]
    #datadel=data.loc[data['Total'] > 2]
    #datadel.to_csv('./result/unbindcleansample/del/'+file+"del.txt",sep='\t',index=False)
    datadeal.to_csv('./result/unbindcleansample/datadeal/'+file+"datadeal.txt",sep='\t',index=False,header=0)
datadel,dataunbind,datadel="","",""

######
if not os.path.exists('./result/unbindcleansample/lastdata/'):
    os.makedirs('./result/unbindcleansample/lastdata/')
for file in sample:
    print ('Start processing file:',file)
    with open("./result/unbindcleansample/datadeal/"+file+"datadeal.txt","r") as pdinf:
        with open("./result/unbindcleansample/lastdata/"+file+"deal.txt","w") as pdout:
            repeat={}
            for lines in pdinf.readlines():
                line=lines.split("\t")
                if line[0] not in repeat.keys():
                    repeat={}
                    linebefore=line
                    repeat[line[0]]=line[0]
                else:
                    start1=list(range(int(linebefore[2]),int(linebefore[2])+len(linebefore[3])))
                    start2=list(range(int(line[2]),int(line[2])+len(line[3])))
                    inter=list(set(start1).intersection(set(start2)))
                    inter1=list(set(start1).difference(set(start2)))
                    inter2=list(set(start2).difference(set(start1)))
                    inter.sort()
                    inter1.sort()
                    inter2.sort()
                    seq1=linebefore[3][inter[0]-int(linebefore[2]):inter[-1]-int(linebefore[2])]
                    seq2=line[3][inter[0]-int(line[2]):inter[-1]-int(line[2])]
                    # data1.loc[:,"class"]="intersection2bind"
                    if len(inter1)!=0 and len(inter2)==0:
                        seqlast=linebefore[3]
                    elif len(inter1)==0 and len(inter2)!=0:
                        seqlast=line[3]
                        linebefore[2]=line[2]
                    elif len(inter1)==0 and len(inter2)==0:
                        seqlast=linebefore[3]
                    else:
                        if inter1[0]<inter2[0]:
                            seqlast=linebefore[3]+line[3][inter2[0]-int(line[2]):]
                        else:
                            seqlast=line[3]+linebefore[3][inter1[0]-int(linebefore[2]):]
                            linebefore[2]=line[2]
                    linebefore[3]=seqlast
                    pdout.write("\t".join(linebefore))  
                        
    # dataunique= pd.read_csv("./result/unbindcleansample/"+file+"unbindlast.txt", names= ["ID","match","start","seq","check"],header=None,sep="\t",index_col=None)                   
    # data= pd.read_csv("./result/unbindcleansample/"+file+"deal.txt", names= ["ID","match","start","seq","check"],header=None,sep="\t",index_col=None)
    # dataunique = pd.concat((dataunique,data),axis=0)
    # dataunique.to_csv('./result/unbindcleansample/'+file+"unbindfinal.txt",sep='\t',index=False,header=0)


                            




