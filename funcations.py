import numpy as np
import bisect
from itertools import permutations

def Complement(seq):
    dic={"G":"C","C":"G","A":"T","T":"A"}
    s = list(seq)
    for i in range(len(seq)):
        s[i]=dic[seq[i]]
    s = "".join(s)
    return s

def Reverse(seq):
    s = list(seq)
    s = reversed(s)
    s = "".join(s)
    return s

def Reverse_Complement(seq):
    seq = Reverse(seq)
    seq = Complement(seq)
    return seq

def Translation_Table(seq):
    dic = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }
    s = ""
    for i in range(0, len(seq)- 2, 3):
        s += dic[seq[i:i+3]]
    return s        

def Match(seq,sub_seq):
    x=-1
    for i in range(len(seq)-len(sub_seq)+1):
        if sub_seq==seq[i:i+len(sub_seq)]:
            x=i
            break
    return x

def Badchars(seq,sub_seq):
    table=np.zeros([4,len(sub_seq)])     
    row=["A","C","G","T"]
    for i in range (4):
        num=-1
        for j in range (len(sub_seq)):
            if row[i]==sub_seq[j]:
                table[i,j]=-1
                num=-1
            else:
                num+=1
                table[i,j]=num
    x=-1
    i=0
    while(i<len(seq)-len(sub_seq)+1):
        if sub_seq==seq[i:i+len(sub_seq)]:
            x=i
            break
        
        else:
            for j in range(len(sub_seq)-1,-1,-1):
                if seq[i+j] != sub_seq[j]:
                    k=row.index(seq[i+j])
                    i+=table[k,j]
                    break
        i=int(i+1)
    return x

def IndexSorted(seq,ln):
    index = []
    for i in range(len(seq)-ln+1):
        index.append((seq[i:i+ln], i))
    index.sort() 
    return index

def query(t,p,index):
    keys = [r[0] for r in index]
    st = bisect.bisect_left(keys,p[:len(keys[0])])
    en = bisect.bisect(keys,p[:len(keys[0])])
    hits = index[st:en] 
    l=[h[1] for h in hits ]
    offsets=[]
    for i in l:
        if t[i:i+len(p)]==p:
            offsets.append(i)
    return offsets

def Suffix_array_construction(T):
    result = []
    dec = {
    '$' : 0,
    'A' : 1,
    'C' : 2,
    'G' : 3,
    'T' : 4
    }
    table=[]
    i=2**0
    n=0
    while True:    
        l=[]
        dec2={}
        if i>1:
            for j in range(len(T)):
                if not(table[n-1][j:j+i] in l):
                    l.append(table[n-1][j:j+i])
            l.sort()
            l
            for j in range(len(l)):
                dec2[tuple(l[j])]=j
        row=[]
        for j in range (len(T)):
            if i==1:
                row.append(dec[T[j]])
            else:
                row.append(dec2[tuple(table[n-1][j:j+i])])
        table.append(row)
        flag=0
        for j in range(len(row)):
            c=0
            c=row.count(j)
            if c>1:
                flag=1
                break
        # print(f'iteration {n}: {row}')
        result.append(f'iteration {n}: {row}')
        if flag==0:
            break
        n+=1
        i=2**n
    return result

def overlap(a, b, min_length=3):
    start=0
    while True:
        start=a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1

def native_overlap(reads, k):
    olap={}
    reads = ''.join(reads.split())
    reads = reads.split(',')
    for a,b in permutations(reads, 2):
        olen=overlap(a, b, k)
        if olen > 0:
            olap[(a, b)]=olen
    return olap

