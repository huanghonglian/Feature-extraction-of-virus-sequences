import argparse
import sys
import re
import os
import time
from concurrent.futures.process import ProcessPoolExecutor
import numpy as np
from scipy.stats import entropy
import pandas as pd
file_name='./public/sample/metadata.txt'
rev_comp = False
REF_FILE = "./public/reference/ref_seq.fasta"
with open(REF_FILE) as fp:
    REF_SEQ = fp.readlines()[1].strip()
    REF_LEN = len(REF_SEQ)
REF_SEQ=REF_SEQ.replace("A","0").replace("C","1").replace("G","2").replace("T","3")
def trans(ch):
    if ch in ['a', 'A']:
        return 0
    if ch in ['c', 'C']:
        return 1
    if ch in ['g', 'G']:
        return 2
    if ch in ['t', 'T']:
        return 3
    return -1

def extract_id(s):
    return s.split('/')[-1].split('.')[0]

class Nuc2aa_converter:
    def __init__(self, filename):
        self.mapping = [0] * 64
        with open(filename) as input_file:
            lines = input_file.readlines()
            self.total = len(lines)
            for i, line in enumerate(lines):
                for codon in line.strip().split(", "):
                    h = trans(codon[0]) * 16 + trans(codon[1]) * 4 + trans(
                        codon[2])
                    self.mapping[h] = i

    def nuc2aa(self, nuc):
        h = nuc[0] * 16 + nuc[1] * 4 + nuc[2]
        return self.mapping[h]

    def nuc2int(self, nuc):
        h = nuc[0] * 16 + nuc[1] * 4 + nuc[2]
        return h
converter = Nuc2aa_converter("../public/reference/DNA_codon_table.txt")

class Features:
    def __init__(self, file_id):
        self.file_id = file_id
        self.max_mut = 1
        self.E6_max_mut = 1
        self.sumcount = np.zeros(REF_LEN,dtype="float")
        self.mutation = np.zeros(REF_LEN,dtype="float")
        self.count = np.zeros((5, REF_LEN),dtype="float")
        self.entropy = []
        for i in range(REF_LEN):
            self.entropy.append([])
        self.group = np.zeros((64, REF_LEN // 3),dtype="float")
        self.slide = np.zeros((64, REF_LEN - 2),dtype="float")

    def statistic(self):
        # print("nuc", file =sys.stderr)
        start_time = time.time()
        self.entropy_nuc = entropy(self.count, base=2)
        self.entropy_nuc=np.where(self.entropy_nuc > 0, self.entropy_nuc, 0)       
        self.mutation = np.divide(self.mutation, self.sumcount)
        self.mutation=np.nan_to_num(self.mutation)
        self.count = np.divide(self.count, self.sumcount) 
        self.count=np.nan_to_num(self.count)
        #print("group", file =sys.stderr)
        self.entropy_group = entropy(self.group, base=2)
        self.entropy_group=np.where(self.entropy_group > 0, self.entropy_group, 0)
        # print("slide", file =sys.stderr)
        self.entropy_slide = entropy(self.slide, base=2)
        self.entropy_slide=np.where(self.entropy_slide > 0, self.entropy_slide, 0)
        self.sumcount  =  self.sumcount 
        end_time = time.time()
        #print("Time2: %d sec" % int(end_time - start_time))

    def out(self):
        out_str = self.file_id
        #print (out_str)
        #out_str=out_str+"\t"+str(sample.loc[out_str,"Host_disease"])
        start_time = time.time()
        count=self.count[0]
        for j in range(1,4):
            count=np.concatenate((count,self.count[j]),axis=0)
        mutation=out_str+"\t"+"\t".join(map(str,self.mutation))
        count=out_str+"\t"+"\t".join(map(str,count))
        entropy_nuc=out_str+"\t"+"\t".join(map(str,self.entropy_nuc))
        entropy_slide=out_str+"\t"+"\t".join(map(str,self.entropy_slide))
        entropy_group=out_str+"\t"+"\t".join(map(str,self.entropy_group))
        sumcount=out_str+"\t"+"\t".join(map(str,self.sumcount))
        end_time = time.time()
        #print("Time3: %d sec" % int(end_time - start_time))
        return mutation,count,entropy_nuc,entropy_slide,entropy_group,sumcount

def parse_fastq(data, f):
    step = 1
    start_time = time.time()
    for k in range(0, len(data), step):
        line = data[k].strip().split("\t")
        seq =line[2].upper()
        seq_trans = [0] * REF_LEN
        seq=seq.replace("A","0").replace("C","1").replace("G","2").replace("T","3")
        seq = re.sub(r'\D', "4",seq)   
        # 1 nucleotide
        for x in range(len(seq)):
            pos=x+int(line[1])-1
            seq_trans[pos] = int(seq[x])
            f.count[seq_trans[pos]][pos] += 1
            f.sumcount[pos]+=1
            if seq_trans[pos]!=4:
                f.entropy[pos].append(seq_trans[pos])
             # 1 mutatation
            if not seq[x] == REF_SEQ[pos] and seq[x]!="4":
                f.mutation[pos] += 1
                if f.mutation[pos] > f.max_mut:
                    f.max_mut = f.mutation[pos]
        # 3 nucleotides (slide)
            if pos < len(seq)-2 and 4 not in seq_trans[pos:pos + 3]:
                kmer = converter.nuc2int(seq_trans[pos:pos + 3])
                f.slide[kmer][pos] += 1
        ##entropy
        # 3 nucleotides (group)
                if len(seq)%3==0 and 4 not in seq_trans[pos:pos + 3]:
                    kmer = converter.nuc2int(seq_trans[pos:pos + 3])
                    f.group[kmer][pos // 3] += 1
    end_time = time.time()
    print("Time1: %d sec" % int(end_time - start_time))

def work(file_name):
    f = Features(extract_id(file_name))
    with open("./result/finaldata/seq/"+file_name+".txt","r") as input_file:
        parse_fastq(input_file.readlines(), f)
        pass
    f.statistic()
    return f.out()
if not os.path.exists('./feature/data/'):
    os.makedirs('./feature/data/')
if not os.path.exists('./feature/data/raw/'):
    os.makedirs('./feature/data/raw/')
if not os.path.exists('./feature/data/each/'):
    os.makedirs('./feature/data/each/')
if __name__ == "__main__":
    sample= pd.read_csv(file_name, header=0,index_col=0,sep="\t")
    file_list = list(sample.index.tolist())
    for l in file_list:
        print ('Start processing file:',l)
        with open("./feature/data/raw/"+l+"mutation.txt", "w") as outfile1:
            with open("./feature/data/raw/"+l+"count.txt", "w") as outfile2:
                with open("./feature/data/raw/"+l+"entropy_nuc.txt", "w") as outfile3:
                    with open("./feature/data/raw/"+l+"entropy_slide.txt", "w") as outfile4:
                        with open("./feature/data/raw/"+l+"entropy_group.txt", "w") as outfile5:
                            with open("./feature/data/each/"+l+"mutation.txt", "w") as outfile6:
                            #with ProcessPoolExecutor(max_workers=40) as executor:
                                res= work(l)     
                                outfile1.write(res[0])
                                outfile1.write("\n")         
                                outfile2.write(res[1])
                                outfile2.write("\n")       
                                outfile3.write(res[2])
                                outfile3.write("\n")             
                                outfile4.write(res[3])
                                outfile4.write("\n")              
                                outfile5.write(res[4])
                                outfile5.write("\n")
                                outfile6.write(res[5])
                                outfile6.write("\n")
