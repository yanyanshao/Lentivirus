# -*- coding:utf-8 -*-
###################################################################
# File Name: star-junction.correct.LTpos.py
# Author: yys
# mail: shayy0919@163.com
# Created Time: 2019年08月11日 星期日 21时06分51秒
###################################################################
#!/usr/bin/python3

import os
import sys
import pdb
from Bio import SeqIO
from math import log
from collections import Counter


def GetLocus(refDBfp):
    genomeDict = {}
    fasta_sequences = SeqIO.parse(refDBfp,'fasta')

    for fasta in fasta_sequences:
        chrid = fasta.id
        chrseq    = str(fasta.seq)
        genomeDict[chrid] = chrseq

    return genomeDict

def ComputeEntropy(seq):
    #setseq = set(seq)
    total = len(seq)
    seqList = list(seq)
    countDict = Counter(seqList)

    entropy = 0
    for key, value in countDict.items():
        p = float(value)/total
        entropy += p * (log(1/p) / log(2))

    return entropy


def ComputeEntropySlideWindow(chrom, genDict, slide):
    seq = genDict[chrom]
    seql = len(seq)

    total = 0
    for i in xrange(seql-slide):
        entroy = ComputeEntropy(seq[i:i+slide])
        #print("%d\t%f" %((i+i+slide)/2,entroy))
        total += entroy
    
    average = float(total)/(seql-slide)
    
    return average


def main():
    #pdb.set_trace()
    juncfp = open("/home/yys/Project/Lilab/Software/Result/star-fusion.multi.map.onebc.pass.txt", "r")
    genomeDict = GetLocus("/home/yys/Project/Lilab/Lentivirus/DataBase/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/ref_genome.fa")
    
    averageDict = {}
    for key, value in genomeDict.items():
        average = ComputeEntropySlideWindow(key, genomeDict, 16)
        averageDict[key] = average

    for line in juncfp.readlines():
        lline = line.rstrip().split()

        splicetype = lline[6]
        if splicetype == '0':
            if lline[0] == "chrLT":
                chrom = lline[3]
                juncseq = genomeDict[lline[3]][int(lline[4])-8:int(lline[4])+8]
            else:
                chrom = lline[0]
                juncseq = genomeDict[lline[0]][int(lline[1])-8:int(lline[1])+8]
        entroy = ComputeEntropy(juncseq)
        print("%s\t%f\t%f" %(line.rstrip(),entroy),averageDict[chrom])
    
    #a = ComputeEntropy("ATGCGTCGTATT")
    #b = ComputeEntropy("TGTCGTCCCTAGT")


if __name__ == '__main__':
    main()
