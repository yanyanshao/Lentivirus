#!/usr/bin/python2.7

import sys
import os
import time
import shutil
import pandas as pd
from parameters import *
from Bio import SeqIO
from sequencefunction import GcContent

t0 = time.time()


'''
    refdb_path[directory]               vcf_path[directory]
            |                                   |
            |                                   |
            |                                   |
    ref_fasta                               poly_vcf
        |                                       |
        |                                       |
 chr1.fasta,chr2.fasta                      chr1,chr2,chr3
'''

''' Function: read whole genome_fasta file
     Input:                           
           @NUll
     Output:                          
           @Null
'''
def GenomeRead(fa):
    fasta_seq = SeqIO.parse(fa, "fasta")
    genome_dict = {}
    for fasta in fasta_seq:
        genome_dict[fasta.id] = [str(fasta.seq).upper(), len(fasta.seq)]    #genome_dict is global variable

    return genome_dict


def main():
    fa_fp = open(cart_fa, "r")
    poly_seq_fp = open(outdir+"/sequence.txt", "w")

    fa_dict = GenomeRead(fa_fp)

    for key, seq in fa_dict.items():
        gc_locus = GcContent(seq[0])
        poly_seq_fp.write(">TA" + "_" + "lentivirus" + "_" + "0"  + "_" + str(seq[1]) + "_ " + "0"  + "_" + "0"  + "_"+ "0"  + "_" + str(gc_locus) + "_" + "1" + "\n")
        poly_seq_fp.write(seq[0] + "\n")
    
    #write parameters used out put file
    parameters_used = open(outdir +'run_summary.txt', 'w')

    #Print parameters used  
    parameters_used.write   (   "##########################################################"+"\n"+
                                "### TRS locus selection & variant masking parameters"+"\n"+
                                "### Start date: " + str(time.strftime("%m/%d/%Y"))+";"+ "\t"+ "Start time: "+ str(time.strftime("%H:%M:%S")) +"\n"+
                                "##########################################################"+"\n")

    t1 = time.time()
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))

    parameters_used.close()

    # Copy parameters file into the timestamp directory
    shutil.copy2('./parameters.py', str(outdir + '/parameters.py'))
               

if __name__ == '__main__':
    main()
