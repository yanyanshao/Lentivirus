#!/usr/bin/python

import sys
import os

def GcContent(sequence):
    sequence = sequence.upper()
    gc_count =  sequence.count('G') + sequence.count('C')
    seq_len = len(sequence)
    gc_cont = (float(gc_count)/seq_len)*100

    return round(gc_cont, 2)
    
def rev_complement(seq):
        seq = seq.upper()
        basecomplement = {'A':'T',
                                          'C':'G',
                                          'G':'C',
                                          'T':'A',
                                          '-':'-',
                                          'N':'N'}
        letters = list(seq)
        letters = [basecomplement[base] for base in letters]
        complement = (''.join(letters))                                                 #gives the complement of the bases in list letters
        return complement[::-1]                                                         #gives the reverse of the compliment

