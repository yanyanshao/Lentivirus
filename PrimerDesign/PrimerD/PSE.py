#!/usr/bin/python

import os
import re
import pdb
import time
import primer3
import pandas as pd
from Bio import SeqIO
import subprocess as sp
import numpy as np
from parameters import *
from sequencefunction import rev_complement, GcContent
from Santalucia_NN_Tm import NN_Tm, complement, mM_monovalent

t0 = time.time()
genome_dict = {}
monovalent_cation_eq    =    mM_monovalent(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs)


''' Function: calculte form hairpin needed tm
     Input: 
           @primer_sequence : DNA primer [string]
     Output: 
           @Tm_hairpin    : form hair needed tm
'''
def hairpin_Tm(primer_sequence, mv_cation=0,primer_conc=0):
    #use primer3 calculate form hairpin tm
    Tm_hairpin =  (primer3.calcHairpin(primer_sequence,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    return ("{0:.2f}".format(round(Tm_hairpin,2)))


''' Function: calculte form homodimer needed tm
     Input: 
           @primer_sequence : DNA primer [string]
     Output: 
           @Tm_homodimer    : form homodimer needed tm
'''
def homodimer_Tm(primer_sequence, mv_cation=0,primer_conc=0):
    Tm_homodimer = (primer3.calcHomodimer(primer_sequence,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    return ("{0:.2f}".format(round(Tm_homodimer,2)))


''' Function: calculte primer tm use normal tm
     Input: 
           @sequence : DNA primer [string]
     Output: 
           @Tm       : primer sequence tm
'''
def Tm(sequence):
    sequence = sequence.upper()
    A = sequence.count('A')
    C = sequence.count('C')
    G = sequence.count('G')
    T = sequence.count('T')
    Tm = (2*(A+T))+ (4*(G+C))   #calculte with normal methods
    return Tm


''' Function: calculate continuous gc stretch
     Input:                           
           @sequence : DNA primer [string]
     Output:                          
           @out_put         : gc stretch num
             @0             : no "G" or "C"
             @1             : have "G" or "C", only one
             @other(2,3...) : have not one ,
'''
def continuous_gc(sequence):
    sequence     = sequence.upper()
    if "G" in sequence or "C" in sequence :
        #instead "C" with "G"
        sequence     = sequence.replace("C","G")
        if "GG" in sequence:
            #calculte "GG" num
            out_put        = len(max(re.compile("(G+G)").findall(sequence)))
        else:
            out_put        = 1
    else:
        out_put        = 0
    return out_put


''' Function: read whole hg19.fasta file
     Input:                           
           @NUll
     Output:                          
           @Null
'''
def GenomeRead():
    #hg19.fasta from parameter.py file(hg19.fasta)
    fasta_seq = SeqIO.parse(genome_fasta, "fasta")
    for fasta in fasta_seq:
        genome_dict[fasta.id] = [str(fasta.seq), len(fasta.seq)]    #genome_dict is global variable
    

''' Function: sequence extraction from the genome
     Input:                           
           @fasta_id :  chromsome name 
           @coordinates_3primeextend : position coordinate needed extract[1035:1200]
     Output:                          
           @select_seq : follow the position sequence
'''
def seq_extraction(fasta_id, coordinates_3primeextend):
    chr_seq     =    fasta_id
    coordinates =    [x.strip() for x in coordinates_3primeextend.split(':')]
    start_coord =    int(coordinates[0])
    stop_coord  =    int(coordinates[1])
    if start_coord < 1:         #if start_coord < 1,it is unnormal position 
        start_n = start_coord
        start_coord = 1
        #extract sequence from genome_dict
        select_seq  =    genome_dict[str(chr_seq)][0][(start_coord-1):(stop_coord)]
        #complement with "N"
        select_seq  =    ('N'*start_n) + select_seq
    elif stop_coord > genome_dict[str(chr_seq)][1]:     #judge end position whether unnormal
        stop_n = stop_coord
        start_coord = int(genome_dict[str(chr_seq)][1])
        select_seq  =    genome_dict[str(chr_seq)][0][(start_coord-1):(stop_coord)]
        #unnormal position comple with base "N"
        select_seq  =    select_seq+('N'*stop_n)
    else:
        select_seq  =    genome_dict[str(chr_seq)][0][(start_coord-1):(stop_coord)]
    return select_seq



''' Function: calculate 3' end mismatch num
     Input:                           
           @seq1 : sequence 1
           @seq2 : sequence 2
                   AGTGCGTGCAGTC
                   AGTGCGTGCAGAT        --return 2
                    
                   AGTGCGTGCAGTC
                   AGTGCGTGCTTTC        --return 0
     Output:                          
           @num_mismatch : 3' end mismatch num
'''
def count_3prime_mismatches(seq1, seq2):
    num_mismatch = 0
    seq1 = (seq1[::-1]).upper()
    seq2 = (seq2[::-1]).upper()
    list = [seq1,seq2]
    for i in xrange(len(min(list, key=len))):
        #if seq base match,break
        if seq1[i] == seq2[i]:
            break
        elif seq1[i] != seq2[i]:
            num_mismatch +=1
    return num_mismatch


''' Function: function to count the start position of mismatches from 3' end of two strands
     Input:                            
           @seq1 : sequence 1
           @seq2 : sequence 2
                   AGTGCGTGCAGTC
                   AGTGCGTGCAGTT        --return 1
                    
                   AGTGCGTGCAGTC
                   AGTGCGTGCATTC        --return 3
     Output:                             
           @mismatches_start_pos3prime : first 3' end mismatch position
'''
def start_pos3prime_mismatches(seq1, seq2):
    mismatches_start_pos3prime = "Nill"
    seq1 = (seq1[::-1]).upper()
    seq2 = (seq2[::-1]).upper()
    list = [seq1,seq2]
    for i in xrange(len(min(list, key=len))):
        if seq1[i] != seq2[i]:
            mismatches_start_pos3prime = i+1
            break
    return mismatches_start_pos3prime


''' Function: function to count the hamming distance between two sequences
     Input:                            
           @seq1 : sequence 1
           @seq2 : sequence 2
                   AGTGCGTGCAGTC
                   AGTGAGTGCAGTT        --return 2
                          
                   AGTGCGTGCAGTC
                   ATTGCGTGCATTA        --return 3
     Output:                             
           @score : two sequence total mismatch num
'''
def HammingDistance(seq1, seq2):
    if len(seq1) == len(seq2):
        score = 0
        seq1 = (seq1).upper()
        seq2 = (seq2).upper()
        for i in xrange(len(seq1)):
            if seq1[i] != seq2[i]:
                score = score+1
    else:
        score = "different_lengths!"
    return score


''' Function: function to calculate the percentile of the mipriming tm s for a particular primer
     Input:                            
           @Tm_list                      : all off-target sequence tm[list]
           @misprime_Tm_percentile_value : percentile[eg:90]
     Output:                             
           @p                            : percentile of the mipriming tm s for a particular primer
'''
def misprime_percentile(Tm_list, misprime_Tm_percentile_value):
    Tm_list_cleaned = [ x for x in Tm_list if isinstance(x, float)]
    if len(Tm_list_cleaned) == 0:
        p = 'No match/3_p_mismatch'
    else:
        #calculte use numpy percentile function
        p = np.percentile(Tm_list_cleaned, misprime_Tm_percentile_value)
        p =(np.float64(p).item())
        p = round(p,2)
    return p


''' Function: function to calculate the max mipriming tm for a particular primer
     Input:                                                                              
           @Tm_list                      : all off-target sequence tm[list]
     Output:                             
           @max_tm                       : max mipriming tm for a particular primer
'''
def max_misprime_Tm(Tm_list):
    Tm_list_cleaned = [ x for x in Tm_list if isinstance(x, float)]
    if len(Tm_list_cleaned) == 0:
        max_Tm = "NA"
    else:
        max_Tm    = max(Tm_list_cleaned)
    return max_Tm


''' Function: function to calculate the min mipriming length for a particular primer
     Input:                                                                              
           @pdict                      : all off-target sequence tm[dict]
             {id:[],id_2:[]...}        : id = ">TA_chr7_....._ATGCGTGCGTGCAGTCGT"
     Output:                                                         
           @len(min(p_list, key=len))  : min mipriming length
'''
def min_primer_length(pdict):
    p_list = []
    for key in pdict.keys():
        primer = key.split("_")[-1]
        p_list.append(primer)
    return len(min(p_list, key=len))

def main():
    three_primer_region = 5
    word_size = 7
    misprime_Tm_percentile_value = 90
    forward_primer = outdir+"UOD_forward_primer.fasta"
    reverse_primer = outdir+"UOD_reverse_primer.fasta"

    parameters_used = open(outdir + 'run_summary.txt', 'a')

    d = dict()
    d2 = dict()

    #read whole hg19.fasta file
    GenomeRead()
    #create two dict for the primer sequence
    with open(forward_primer, "r") as f:
        for line in f.readlines():
            if  line.startswith(">"):
                d[line[1:-1]] = ['No match',0, 'No match', 0]
                d2[line[1:-1]] = []

    #get the min primer length in all potential primer sequence
    min_primer_f_len = min_primer_length(d)
    min_primer_len = min_primer_f_len
    
    #if set 3' length high then min primer length,then raise error
    if three_primer_region > min_primer_len:
        raise ValueError('3 primer region should not greater than the mininum primer length! Given three_primer_region = ', three_primer_region, 'min_primer_len = ', min_primer_len)

    f0 = open(os.devnull, "w")
    #blast to find new similar matchs for forward primer sequence
    p1 = sp.Popen(["blastn","-task","blastn","-db","%szmv2all" %refDB_path,"-query","%s" %(forward_primer),"-evalue","%s" %mp_e_value,"-word_size","%s" %word_size,"-gapopen","%s" %mp_gapopen,"-gapextend","%s" %mp_gapextend,"-reward","%s" %mp_reward,"-penalty","%s" %mp_penalty,"-dust","no","-perc_identity","%s" %mp_perc_identity,"-max_target_seqs", "%s" %mp_max_target_seqs,"-max_hsps", "%s" %mp_max_hsps,"-outfmt","10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq qseq", "-num_threads","%s" %mp_num_threads], stdout=sp.PIPE,stderr=f0)
    output, error = p1.communicate()

     ### Print parameters used
    parameters_used.write   (  "##########################################################"+"\n"+
                    "### Summary of PSE parameters"+"\n"+
                    "### Start date: " + str(time.strftime("%m/%d/%Y"))+";"+ "\t"+ "Start time: "+ str(time.strftime("%H:%M:%S"))+"\n"+
                    "##########################################################"+"\n"+
                    "## blastn parameters for mis-prime alignment"+"\n"+
                    "Query_file_F_primers       =   "+str(forward_primer)+"\n"+
                    "Query_file_R_primers       =   "+str(reverse_primer)+"\n"+
                    "e_value                    =   "+str(mp_e_value)+"\n"+
                    "Word_size_F_primers        =   "+str(word_size)+"\n"
                            )

    #manage the blast result
    list_lines = []
    for line in output.split('\n')[:-1]:
        line = line.rstrip(' ').split(',')
        query_len   =    len(line[0].split("_")[-1])
        qstart        =    int(line[6])
        qend        =    int(line[7])
        sstart        =    float(line[8])
        send        =    float(line[9])
        pident        =    float(line[2])
        match_length    =    float(line[3])

        #exact match sequence
        if pident == 100 and query_len == match_length:
            pass
        #off-target sequence
        else:
            three_primer_mismatchs = query_len - qend   #3' end mismatch num
            five_primer_mismatchs = qstart - 1          #5' start mismatch num
            fasta_id = line[1]                  #chromsome
            Match_seq = line[12]                # match seq
            Query_seq = line[13]                #query seq
            match_gaps = Match_seq.count('-')   #match gaps
            query_gaps = Query_seq.count('-')   #query gaps
            if send-sstart > 0:
                #strand "+"
                match_direction = 'right'
                actual_match_start = int(sstart)
                actual_match_end = int(send)
                #recalculate the position of sequence
                actual_match_end_3prime_extended = actual_match_end + int(three_primer_mismatchs)
                actual_match_start_5prime_extended = actual_match_start - int(five_primer_mismatchs)
                actual_match_start_5prime_extended = actual_match_start_5prime_extended + query_gaps
                actual_match_start_5prime_extended = actual_match_start_5prime_extended - match_gaps
                if int(three_primer_mismatchs) == 0 and int(five_primer_mismatchs) == 0 and match_gaps == 0 and query_gaps == 0:
                    match_extend_3_5 = Match_seq
                else:
                    #form new position coordinate
                    coordinates_3_5_extended = str(actual_match_start_5prime_extended) + ":" + str(actual_match_end_3prime_extended)
                    #extract new sequence
                    match_extend_3_5 = seq_extraction(fasta_id, coordinates_3_5_extended)
            else:
                #strand "-"
                match_direction = 'left'
                actual_match_start = int(send)
                actual_match_end = int(sstart)
                actual_match_end_3prime_extended = actual_match_start - int(three_primer_mismatchs)
                actual_match_start_5prime_extended = actual_match_end + int(five_primer_mismatchs)
                actual_match_start_5prime_extended = actual_match_start_5prime_extended + match_gaps
                actual_match_start_5prime_extended = actual_match_start_5prime_extended - query_gaps
                if int(three_primer_mismatchs) == 0 and int(five_primer_mismatchs) == 0 and match_gaps == 0 and query_gaps == 0:
                    match_extend_3_5 = rev_complement(Match_seq)
                else:
                    coordinates_3_5_extended = str(actual_match_end_3prime_extended) + ":" + str(actual_match_start_5prime_extended)
                    match_extend_3_5 = seq_extraction(fasta_id, coordinates_3_5_extended)
            Primer = line[0].split("_")[-1]
            line.extend([match_direction])
            line.extend([three_primer_mismatchs])
            line.extend([five_primer_mismatchs])
            #calculte tm use normal function
            Primer_tm = Tm(Primer)
            line.append(Primer_tm)
            match_melt_temp = Tm(Match_seq)
            line.append(match_melt_temp)
            line.append(Primer_tm-match_melt_temp)
            fasta_id = line[1]
            qstart = int(line[6]) - 1
            qend = int(line[7])
            Match_seq_blast =  line[12]
            line.append(match_extend_3_5)

            #calculte 3' end primer mismatch num and first mismatch position
            if match_direction == 'right':
                actual_3primer_mismatchs = count_3prime_mismatches(Primer, match_extend_3_5)
                mismatch_3primer_start_pos = start_pos3prime_mismatches(Primer, match_extend_3_5)
            elif match_direction == 'left':
                actual_3primer_mismatchs = count_3prime_mismatches(Primer, rev_complement(match_extend_3_5))
                mismatch_3primer_start_pos = start_pos3prime_mismatches(Primer, rev_complement(match_extend_3_5))
            elif "-" in Match_seq_blast:
                actual_3primer_mismatchs = "gap_in_match"
                mismatch_3primer_start_pos = "gap_in_match"
            else:
                actual_3primer_mismatchs = "check_your_sequences!"
                mismatch_3primer_start_pos = "check_your_sequences!"

            line.append(actual_3primer_mismatchs)
            line.append(mismatch_3primer_start_pos)
            #use NN-tm methods calculate tm
            Primer_NN_Tm = NN_Tm(seq=Primer, compl_seq=complement(Primer), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            if match_direction == "right":
                complementary_match    =    complement(match_extend_3_5)
            else:
                complementary_match    =    match_extend_3_5[::-1]

            #Primer and its complementry sequence Tm
            Gap_adjusted_end_filling_Tm = NN_Tm(seq=Primer, compl_seq=complementary_match, primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)

            # Query sequence ant Match sequence Tm 
            Local_alignment_Tm = NN_Tm(seq=Query_seq, compl_seq=complement(Match_seq), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            line.append(Primer_NN_Tm)
            line.append(Local_alignment_Tm)
            line.append(Gap_adjusted_end_filling_Tm)
            list_lines.append(line)

            '''function : fill the result dict d&d2
                      d : {>TA_chr7_1000_..._AGTGCGTGCAGTGGCCGT:["3prime_mismatch",3' end have mismatch primer num,NA,total num],}   
                      d2: {>TA_chr7_1000_..._AGTGCGTGCAGTGGCCGT:[every off-target tm]}
            '''
            if actual_3primer_mismatchs > 0 and(d[line[0]][0] == "3prime_mismatch" or d[line[0]][0] == "No match"):
                d[line[0]][0] =    "3prime_mismatch"
            if actual_3primer_mismatchs == 0:
                if (d[line[0]][0]=="No match" )or (d[line[0]][0]=="3prime_mismatch") or (float(Gap_adjusted_end_filling_Tm) > float(d[line[0]][0])) :
                    d[line[0]][0] =    Gap_adjusted_end_filling_Tm
                if float(Gap_adjusted_end_filling_Tm) > float(d[line[0]][0]) :
                    d[line[0]][3] = 1
                    if mismatch_3primer_start_pos > three_primer_region:
                        d[line[0]][1]    += 1
                if  float(Gap_adjusted_end_filling_Tm) == float(d[line[0]][0]) :
                    d[line[0]][3] += 1
                    if mismatch_3primer_start_pos > three_primer_region:
                        d[line[0]][1]    += 1
            d2[line[0]].append(float(Gap_adjusted_end_filling_Tm))

    for key, value in d2.iteritems():
        #get off-target list percentile value, put d dict
        d[key][2]    =    misprime_percentile(d2[key], misprime_Tm_percentile_value)
        #get the max off-target tm
        maximum_misprime_Tm    =    max_misprime_Tm(d2[key])
        if maximum_misprime_Tm != 'NA' and d[key][0] != "3prime_mismatch":
            d[key][0]    =    maximum_misprime_Tm

    headers = ["Primer","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","e-value","bitscore","Match_seq","Query_seq","match_direction","3prime_overhang","5prime_overhang","PrimerTm","MatchTm","PrimerTm-MatchTm", "3&5_prime_match_extend",  "Actual_3prime_mismatches","First3Prime_mismatch" ,"Primer_NN_Tm","Local_alignment_Tm", "Gap_adjusted_end_filling_Tm"]        
    df = pd.DataFrame(list_lines, columns=headers)


    output_headers = ["Primer","qstart","qend","sstart","send","Query_seq","Match_seq","match_direction","3prime_overhang","5prime_overhang","3&5_prime_match_extend","PrimerTm","MatchTm","PrimerTm-MatchTm","First3Prime_mismatch","Primer_NN_Tm","Local_alignment_Tm", "Gap_adjusted_end_filling_Tm"]
    df.to_csv(outdir+ 'PSE_after_blast_forward.csv', columns = output_headers)
    
    primer_f_dict = dict()
    for key,value1 in d.iteritems():
        p = key.split("_")[-1]
        value = value1[0]
        three_prime_mismatch_alignments = value1[3] - value1[1]
        max_tm_alignments = value1[3]       #the largest aligmnet number
        three_primer_mismatch_alignments  = str(three_prime_mismatch_alignments) + "-/-" + str(max_tm_alignments)
        misprime_Tm_percentile = value1[2]
        key_Tm = NN_Tm(seq=p, compl_seq=complement(p), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
        if value != "No match" and value != "3prime_mismatch":
            Tm_difference = str(float(key_Tm)-float(value))
        else:
            Tm_difference = "-"
        TmHairpin    =    hairpin_Tm(p, monovalent_cation_eq, primer_conc)
        TmHomodimer    =    homodimer_Tm(p, monovalent_cation_eq, primer_conc)
        primer_f_dict[key.upper()]    =    [key_Tm, value, Tm_difference, three_primer_mismatch_alignments, misprime_Tm_percentile,TmHairpin,TmHomodimer]

    f.close()


    #######################################################
    #########           reverse primer          ###########
    #######################################################

    d = dict()
    d2 = dict()

    with open(reverse_primer, "r") as r:
        for line in r.readlines():
            if  line.startswith(">"):
                d[line[1:-1]] = ['No match',0, 'No match', 0]
                d2[line[1:-1]] = []

    min_primer_r_len = min_primer_length(d)
    if min_primer_r_len < min_primer_f_len:
        min_primer_len = min_primer_r_len

    word_size = 7
    if three_primer_region > min_primer_len:
        raise ValueError('3 prime region should not be greater than the minimum primer length! Given three_prime_region = ', three_primer_region, 'min_pimer_len = ', min_primer_r_len)    

    f0 = open(os.devnull, 'w')
    p2 = sp.Popen(["blastn","-task","blastn","-db","%szmv2all" %refDB_path,"-query","%s" %(reverse_primer),"-evalue","%s" %mp_e_value,"-word_size","%s" %word_size,"-gapopen","%s" %mp_gapopen,"-gapextend","%s" %mp_gapextend,"-reward","%s" %mp_reward,"-penalty","%s" %mp_penalty,"-dust","no","-perc_identity","%s" %mp_perc_identity,"-max_target_seqs", "%s" %mp_max_target_seqs,"-max_hsps", "%s" %mp_max_hsps,"-outfmt","10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq qseq", "-num_threads","%s" %mp_num_threads], stdout=sp.PIPE,stderr=f0)
    output, error = p2.communicate()

    parameters_used.write(
                    "Word_size_R_primers        =   "+str(word_size)+"\n"+
                    "Gapopen                    =   "+str(mp_gapopen)+"\n"+
                    "Gapextend                  =   "+str(mp_gapextend)+"\n"+
                    "Reward                     =   "+str(mp_reward)+"\n"+
                    "Penalty                    =   "+str(mp_penalty)+"\n"
                    "%_identity                 =   "+str(mp_perc_identity)+"\n"+
                    "Max_target_seqs            =   "+str(mp_max_target_seqs)+"\n"+
                    "Max_hsps                   =   "+str(mp_max_hsps)+"\n"+
                    "Number_threads             =   "+str(mp_num_threads)+"\n"+"\n"+
                    "## PSE parameters"+"\n"+
                    "3_prime_region             =   "+str(three_primer_region)+"\n"+
                    "Misprime_Tm_percentile     =   "+str(misprime_Tm_percentile_value )+"\n"+
                    "##########################################################"+"\n"
                        )

    list_lines = []
    for line in output.split('\n')[:-1]:
        line = line.rstrip(' ').split(',')
        query_len   =    len(line[0].split("_")[-1])
        qstart        =    int(line[6])
        qend        =    int(line[7])
        sstart        =    float(line[8])
        send        =    float(line[9])
        pident        =    float(line[2])
        match_length    =    float(line[3])

        if pident == 100 and query_len == match_length:
            pass
        else:
            three_primer_mismatchs = query_len - qend
            five_primer_mismatchs = qstart - 1
            fasta_id = line[1]
            Match_seq = line[12]
            Query_seq = line[13]
            match_gaps = Match_seq.count('-')
            query_gaps = Query_seq.count('-')
            if send-sstart > 0:
                match_direction = 'right'
                actual_match_start = int(sstart)
                actual_match_end = int(send)
                actual_match_end_3prime_extended = actual_match_end + int(three_primer_mismatchs)
                actual_match_start_5prime_extended = actual_match_start - int(five_primer_mismatchs)
                actual_match_start_5prime_extended = actual_match_start_5prime_extended - match_gaps
                actual_match_start_5prime_extended = actual_match_start_5prime_extended + query_gaps
                if int(three_primer_mismatchs) == 0 and int(five_primer_mismatchs) == 0 and match_gaps == 0 and query_gaps == 0:
                    match_extend_3_5 = Match_seq
                else:
                    coordinates_3_5_extended = str(actual_match_start_5prime_extended) + ":" + str(actual_match_end_3prime_extended)
                    match_extend_3_5 = seq_extraction(fasta_id, coordinates_3_5_extended)
            else:
                match_direction = 'left'
                actual_match_start = int(send)
                actual_match_end = int(sstart)
                actual_match_end_3prime_extended = actual_match_start - int(three_primer_mismatchs)
                actual_match_start_5prime_extended = actual_match_end + int(five_primer_mismatchs)
                actual_match_start_5prime_extended = actual_match_start_5prime_extended + match_gaps
                actual_match_start_5prime_extended = actual_match_start_5prime_extended - query_gaps
                if int(three_primer_mismatchs) == 0 and int(five_primer_mismatchs) == 0 and match_gaps == 0 and query_gaps == 0:
                    match_extend_3_5 = rev_complement(Match_seq)
                else:
                    coordinates_3_5_extended = str(actual_match_end_3prime_extended) + ":" + str(actual_match_start_5prime_extended)
                    match_extend_3_5 = seq_extraction(fasta_id, coordinates_3_5_extended)
            
            Primer = line[0].split("_")[-1]
            line.extend([match_direction])
            line.extend([three_primer_mismatchs])
            line.extend([five_primer_mismatchs])
            Primer_tm = Tm(Primer)
            line.append(Primer_tm)
            match_melt_temp = Tm(Match_seq)
            line.append(match_melt_temp)
            line.append(Primer_tm-match_melt_temp)
            fasta_id = line[1]
            qstart = int(line[6]) - 1
            qend = int(line[7])
            Match_seq_blast =  line[12]
            line.append(match_extend_3_5)

            if match_direction == 'right':
                actual_3primer_mismatchs = count_3prime_mismatches(Primer, match_extend_3_5)
                mismatch_3primer_start_pos = start_pos3prime_mismatches(Primer, match_extend_3_5)
            elif match_direction == 'left':
                actual_3primer_mismatchs = count_3prime_mismatches(Primer, rev_complement(match_extend_3_5))
                mismatch_3primer_start_pos = start_pos3prime_mismatches(Primer, rev_complement(match_extend_3_5))
            elif "-" in Match_seq_blast:
                actual_3primer_mismatchs = "gap_in_match"
                mismatch_3primer_start_pos = "gap_in_match"
            else:
                actual_3primer_mismatchs = "check_your_sequences!"
                mismatch_3primer_start_pos = "check_your_sequences!"

            line.append(actual_3primer_mismatchs)
            line.append(mismatch_3primer_start_pos)
            Primer_NN_Tm = NN_Tm(seq=Primer, compl_seq=complement(Primer), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            if match_direction == "right":
                complementary_match    =    complement(match_extend_3_5)
            else:
                complementary_match    =    match_extend_3_5[::-1]

            #Primer and its complementry sequence Tm
            Gap_adjusted_end_filling_Tm = NN_Tm(seq=Primer, compl_seq=complementary_match, primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            # Query sequence ant Match sequence Tm 
            Local_alignment_Tm = NN_Tm(seq=Query_seq, compl_seq=complement(Match_seq), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
            line.append(Primer_NN_Tm)
            line.append(Local_alignment_Tm)
            line.append(Gap_adjusted_end_filling_Tm)
            list_lines.append(line)

            if actual_3primer_mismatchs > 0 and(d[line[0]][0] == "3prime_mismatch" or d[line[0]][0] == "No match"):
                d[line[0]][0] =    "3prime_mismatch"
            if actual_3primer_mismatchs == 0:
                if (d[line[0]][0]=="No match" )or (d[line[0]][0]=="3prime_mismatch") or (float(Gap_adjusted_end_filling_Tm) > float(d[line[0]][0])) :
                    d[line[0]][0] =    Gap_adjusted_end_filling_Tm
                if float(Gap_adjusted_end_filling_Tm) > float(d[line[0]][0]) :
                    d[line[0]][3] = 1
                    if mismatch_3primer_start_pos > three_primer_region:
                        d[line[0]][1]    += 1
                if  float(Gap_adjusted_end_filling_Tm) == float(d[line[0]][0]) :
                    d[line[0]][3] += 1
                    if mismatch_3primer_start_pos > three_primer_region:
                        d[line[0]][1]    += 1
            d2[line[0]].append(float(Gap_adjusted_end_filling_Tm))

    for key, value in d2.iteritems():
        d[key][2]    =    misprime_percentile(d2[key], misprime_Tm_percentile_value)
        maximum_misprime_Tm    =    max_misprime_Tm(d2[key])
        if maximum_misprime_Tm != 'NA' and d[key][0] != "3prime_mismatch":
            d[key][0]    =    maximum_misprime_Tm    

    headers = ["Primer","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","e-value","bitscore","Match_seq","Query_seq","match_direction","3prime_overhang","5prime_overhang","PrimerTm","MatchTm","PrimerTm-MatchTm", "3&5_prime_match_extend",  "Actual_3prime_mismatches","First3Prime_mismatch" ,"Primer_NN_Tm","Local_alignment_Tm", "Gap_adjusted_end_filling_Tm"]
    df = pd.DataFrame(list_lines, columns=headers)


    output_headers = ["Primer","qstart","qend","sstart","send","Query_seq","Match_seq","match_direction","3prime_overhang","5prime_overhang","3&5_prime_match_extend","PrimerTm","MatchTm","PrimerTm-MatchTm","First3Prime_mismatch","Primer_NN_Tm","Local_alignment_Tm", "Gap_adjusted_end_filling_Tm"]
    df.to_csv(outdir+ 'PSE_after_blast_reverse.csv', columns = output_headers)

    primer_r_dict = dict()
    for key,value1 in d.iteritems():
        p = key.split("_")[-1]
        value = value1[0]
        three_prime_mismatch_alignments = value1[3] - value1[1]
        max_tm_alignments = value1[3]       #the largest aligmnet number
        three_primer_mismatch_alignments  = str(three_prime_mismatch_alignments) + "-/-" + str(max_tm_alignments)
        misprime_Tm_percentile = value1[2]
        key_Tm = NN_Tm(seq=p, compl_seq=complement(p), primer_conc=primer_conc, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True)
        if value != "No match" and value != "3prime_mismatch":
            Tm_difference = str(float(key_Tm)-float(value))
        else:
            Tm_difference = "-"
        TmHairpin    =    hairpin_Tm(p, monovalent_cation_eq, primer_conc)
        TmHomodimer    =    homodimer_Tm(p, monovalent_cation_eq, primer_conc)
        primer_r_dict[key.upper()]    =    [key_Tm, value, Tm_difference, three_primer_mismatch_alignments, misprime_Tm_percentile,TmHairpin,TmHomodimer]

    r.close()

    ### Combine two dictionaries (primer_f_dict & primer_r_dict) 
    pooled_primer_f_r_dict    =    dict()
    pooled_primer_f_r_dict.update(primer_f_dict)
    pooled_primer_f_r_dict.update(primer_r_dict)  

    #####################################################################
    ###########          coordinate identification          #############
    #####################################################################

    fasta_input_FR = outdir + "UOD_final_all_primer_fasta.fasta"
    f = open(outdir + "PSE_final_result.csv", "w")
    f.write("id" + ',' + "Primer"+','+ "chrom" + ',' + "start"+','+"stop"+','+"Strand"+','+"Primer_Tm"+','+"Max_misprime_Tm"+','+"Tm_difference"+','+'Misprime_Tm_'+str(misprime_Tm_percentile_value)+'th_percentile'+','+"Primer_GC"+','+"Continuous_GC"+','+"3'_region_mismatches"+','+"Hairpin_Tm"+','+"Homodimer_Tm"+'\n')

    word_size = min_primer_len

    f0 = open(os.devnull, 'w')
    p3 = sp.Popen(["blastn","-db","%szmv2all" %refDB_path,"-query","%s" %(fasta_input_FR),"-evalue","0.1","-word_size","%s" %word_size,"-gapopen","0","-gapextend","2","-reward","1","-penalty","-3","-dust","no","-perc_identity","100","-max_target_seqs", "13","-outfmt","10 qseqid length sstart send", "-num_threads","%s" %mp_num_threads],stdout=sp.PIPE,stderr=f0)
    exact_match_output_pooled_primers, error = p3.communicate()

    #after mis-blast exact match primer sequence
    for exact_match_pooled_output_line in exact_match_output_pooled_primers.split("\n")[:-1]:
        exact_match_pooled_output_line = exact_match_pooled_output_line.strip(' ').split(",")
        primer_pooled = exact_match_pooled_output_line[0].split("_")[-1]
        primer_chrom = exact_match_pooled_output_line[0].split("_")[2]  #chromsome
        primer_gc = GcContent(primer_pooled)    #gc content
        primer_continuous_gc = continuous_gc(primer_pooled)
        query_len_pooled = int(len(primer_pooled))
        match_len_pooled = int(exact_match_pooled_output_line[1])
        chrom = exact_match_pooled_output_line[0].split("_")[2]
        sstart = exact_match_pooled_output_line[0].split("_")[3]
        ssend = int(sstart) + query_len_pooled
        strand = exact_match_pooled_output_line[0].split("_")[-2]
        if query_len_pooled == match_len_pooled:
            Primer_Tm_p = exact_match_pooled_output_line[0].split("_")[-3]  #primer tm
            Max_misprime_Tm_p = pooled_primer_f_r_dict[exact_match_pooled_output_line[0].upper()][1]    #get max mismatch primer tm
            Tm_difference_p = pooled_primer_f_r_dict[exact_match_pooled_output_line[0].upper()][2]      #tm difference
            three_prime_region_mismatches = pooled_primer_f_r_dict[exact_match_pooled_output_line[0].upper()][3]    #three primer region mismatchs
            misprime_Tm_percentile = pooled_primer_f_r_dict[exact_match_pooled_output_line[0].upper()][4]
            Hairpin_Tm            = pooled_primer_f_r_dict[exact_match_pooled_output_line[0].upper()][5]    #hairpin tm
            Homodimer_Tm        = pooled_primer_f_r_dict[exact_match_pooled_output_line[0].upper()][6]  #homodimer tm
            f.write(exact_match_pooled_output_line[0] + ',' + str(primer_pooled) +','+ primer_chrom + ',' + str(sstart)+','+str(ssend)+','+str(strand)+','+str(Primer_Tm_p)+','+str(Max_misprime_Tm_p)+','+str(Tm_difference_p)+','+str(misprime_Tm_percentile)+','+str(primer_gc)+','+str(primer_continuous_gc)+','+str(three_prime_region_mismatches)+','+str(Hairpin_Tm)+','+str(Homodimer_Tm)+'\n')

    ############################################################
    #Time to run the code: end timer
    ############################################################
    t1 = time.time()
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))
    parameters_used.write(
                                "### PSE run duration : " + str(total) + " seconds"+'\n'
                                "##########################################################"+"\n"+
                                "\n"+"\n"
                    )

    parameters_used.close()

if __name__ == '__main__':
    main()
