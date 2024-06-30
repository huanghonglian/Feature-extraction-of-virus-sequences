import argparse
from scipy import stats
from concurrent.futures.process import ProcessPoolExecutor
import numpy as np
import glob
import os
import scipy.stats as test
from scipy.stats.mstats import kruskal
import pandas as pd
###数据进行转置
# #一行一样读数据
import glob
file=glob.glob("./feature/data/each/*.txt")
REF_FILE = "./public/reference/ref_seq.fasta"
if not os.path.exists('./feature/data/file/'):
    os.makedirs('./feature/data/file/')
if not os.path.exists('./feature/data/sampledata/'):
    os.makedirs('./feature/data/sampledata/')
####设置测序深度
print('calculating feature matrices for different depths...')
label=[10,30,50,75,100,150,200,250,300,350,400,450,500]
for kk in label:
    label=["mutation","entropy_nuc"]
    #list1=set(list(range(1,171824)))
    for l in label:
        for v in file:	
            site={}
            with open("./feature/data/file/"+str(kk)+v.replace("./feature/data/each/","").replace("mutation.txt","")+l+".txt","w") as inf3:
                with open(v,"r") as inf2:
                    control=inf2.readlines()[0].strip("\n").split("\t")
                    site[control[0]]=control[1:]
                site=pd.DataFrame(site)
                site[control[0]]=site[control[0]].values.astype('float32')
                nub=site[site[control[0]]<kk]
                #print (nub.index)
                v=v.replace("./feature/data/each/","").replace("mutation.txt","")
                #print(v)
                data={}
                with open('./feature/data/raw/'+v+l+".txt","r") as inf1:
                    k=inf1.readlines()[0].strip("\n").split("\t")
                    data[k[0]]=k[1:]
                data=pd.DataFrame(data)
                data.iloc[nub.index,:]=0
                inf3.write(control[0]+"\t")
                inf3.write("\t".join(map(str,data[control[0]].tolist())))
###mutation、entropy_nuc样本文件
print('merging all feature data for samples...')
label=[10,30,50,75,100,150,200,250,300,350,400,450,500]
sample_name= pd.read_csv('./public/sample/metadata.txt', header=0,index_col=None)
sample_name=sample_name["sample"].tolist()
for kk in label:
    sample1= pd.read_csv('./public/sample/sampleinfo.txt', header=0,index_col=None,sep="\t")
    label=["mutation","entropy_nuc"]
    a=[0,1,2,3,4,5]
    for l in label:
        data={}
        for k1 in sample_name:
            host_disease=sample1[sample1["Run"]==k1]
            with open('./feature/data/file/'+str(kk)+k1+l+".txt","r") as pdout:
                value=pdout.readlines()
                for k in value:
                    data[k1]=k.strip().split("\t")
        data=pd.DataFrame(data)
        data=data.T
        data.to_csv('./feature/data/sampledata/'+str(kk)+"_"+l+"_feature.txt",sep='\t',header=None,index=None)
####"mutation"+"entropy_nuc"
sample1= pd.read_csv('./public/sample/sampleinfo.txt', header=0,index_col=None,sep="\t")
label=[10,30,50,75,100,150,200,250,300,350,400,450,500]
for kk in label:
    i=''
    data={}
    with open('./feature/data/sampledata/'+str(kk)+"_"+str(i)+"mutation_feature.txt","r") as pdout:
        value=pdout.readlines()
        for k in value:
            data[k.strip().split("\t")[0]]=k.strip().split("\t")
    data1=pd.DataFrame(data)
    data={}
    with open('./feature/data/sampledata/'+str(kk)+"_"+str(i)+"entropy_nuc_feature.txt","r") as pdout:
        value=pdout.readlines()
        for k in value:
            data[k.strip().split("\t")[0]]=k.strip().split("\t")[1:]
    data=pd.DataFrame(data)
    data=pd.concat([data1,data],axis=0)
    data=data.T
    #print (data.shape)
    data.to_csv('./feature/data/sampledata/'+str(kk)+"_"+str(i)+"me_feature.txt",sep='\t',header=None,index=None)




