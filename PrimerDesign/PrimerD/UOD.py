#!/usr/bin/python2.7

import os
import sys
import pdb
import math
import glob
import time
import primer3
from parameters import *
import subprocess as sp
from operator import itemgetter
from Santalucia_NN_Tm import NN_Tm, complement, mM_monovalent
from sequencefunction import GcContent, rev_complement

t0 = time.time()

''' Function: check if a string is composition of DNA primer[ATGC]
     Input: 
           @primer : DNA primer [string]
     Output: 
           @True   : the DNA primer is only composition ATGC
           @ False : the DNa primer is composition NATGC
'''
def CheckCompostion(primer):
    allow_chars = set("ATGC")

    if set(primer.upper()) <= allow_chars:
        return True
    else:
        return False


'''Function: check the repeat primer 
    Input: 
          @primer           : the primer primer
          @flag             : whether to check have repeat seq
          @threshold        : repeat num
    Output: 
          @filter_condition : whether the primer meets the standard
'''
def NucleotideRepeatFilter(primer, NucleotideRepeatFilter_flag, NucleotideRepeatFilter_threshold):
    if NucleotideRepeatFilter_flag:
        #whether to filter repeat seq
        filter_condition = True
        set_seq_bases = set(primer)     #uniqe primer seq

        #check single base repeat
        for base in set_seq_bases:
            if (base * NucleotideRepeatFilter_threshold) in primer:
                #have repeat seq in primer"AAAA"
                filter_condition = False
                break

        #check double base repeat
        seq_len = len(primer)
        for i in range(0, (seq_len-1)):
            if (primer[i:i+2] * NucleotideRepeatFilter_threshold) in primer:
                filter_condition = False
                break
    else:
        filter_condition = True

    return filter_condition


'''Functionc : calculate hainpin tm 
    Input: 
          @primer           : the primer primer
          @mv_cation        : whether to check have repeat seq
          @primer_conc      : repeat num
    Output: 
          @tm_hairpin : form hairpin needed tm
'''
def HairPinTm(primer, mv_cation=0, primer_conc=0):
    #
    tm_hairpin = (primer3.calcHairpin(primer,mv_conc=mv_cation, \
                  dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    
    return ("{0:.2f}".format((round(tm_hairpin, 2))))


'''Functionc : calculate hainpin tm 
    Input:                      
          @primer           : the primer primer
          @mv_cation        : whether to check have repeat seq
          @primer_conc      : repeat num
    Output: 
          @tm_homodimer : form homodimer needed tm
'''
def HomodimerTm(primer, mv_cation=0, primer_conc=0):
    tm_homodimer = (primer3.calcHomodimer(primer,mv_conc=mv_cation, \
                    dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm

    return ("{0:.2f}".format((round(tm_homodimer, 2))))


'''Functionc : check whether primer end with "G" or "C" 
    Input:                      
          @primer           : the primer primer
          @CheckATends      : whether to check have repeat seq
    Output: 
          @True : the primer sequence last one base are "G" or "C"
'''
def CheckATends(primer, CheckATends_flag):
    if CheckATends_flag:
        primer = primer.upper()
        #instead the last one "C" base of "G"
        primer = primer[:-1] + primer[-1].replace("C","G")
        if primer[-1] == 'G':
            return True
        else:
            return False
    else:
        return True

'''Functionc : check whether primer 3' end has "GC" clamp
    Input:                      
          @primer               : the primer primer
          @CheckGCclamp_flag    : whether to check GC clamp
    Output: 
          @True : the primer sequence last five base have >3 base "G" and "C"
'''
def CheckGCclamp(primer, CheckGCclamp_flag):
    if CheckGCclamp_flag:
        primer = primer.upper()
        clamp = primer[-5:].replace("C","G")
        if clamp.count("G") >= 3:
            return False
        else:
            return True
    else:
        return True


'''Function: check primer whether meet conditions,such as gccontent, tm...
    Input:
          @num                              : primer belong to region num
          @chrom                            : primer chromsome
          @start                            : bed file whole region start
          @ploc                             : start from the whole region start
          @strand                           : "+" or "-"
          @primer_dict                      : the container for every primer 
          @minTm&maxTm                      : parameters from parameters.py about primer tm ranges(60-70)
          @GC_range_min&GC_range_max        : parameters from parameters.py about primer GC content ranges(40-60)
          @CheckGCclamp&CheckATends         : whether to check GC clamp and whether to check AT ends (1 or 0)
          @NucleotideRepeatFilter_flag      : whether to check repeat sequence(1 or 0)[AAAA]
          @NucleotideRepeatFilter_threshold : repeat nucleotide num (3: AAA or TTT or GGG or CCC)
          @self_Tmdiff                      : parameter from parameter.py about tm self_Tmdiff 
          @monovalent_cation_eq:            : from the function mM_monovalent
    Output:
          @primer_dict                      : contain the primer meet conditon
              {"ATGCGTGCAGTCAGCGTGCA":[1,chr7,1000,32,20,1,57,"+"],...}
'''
def PrimerFilter(num, primer, chrom, start, ploc, strand, primer_dict, minTm, maxTm, GC_range_min, GC_range_max, CheckATends_flag, \
                 CheckGCclamp_flag, NucleotideRepeatFilter_flag, NucleotideRepeatFilter_threshold, self_Tmdiff, monovalent_cation_eq):
    #check primer sequence compositon whether contain N or n
    if CheckCompostion(primer):
        #get primer gc content
        primer_gc = GcContent(primer)
        #to determine whether gc content within a given ranges
        if primer_gc >= GC_range_min and primer_gc <= GC_range_max:
            #check AT ends and GC clamp at 3' end
            if CheckATends(primer, CheckATends_flag) and CheckGCclamp(primer, CheckGCclamp_flag):
                #check repeat requence in primer
                if NucleotideRepeatFilter(primer, NucleotideRepeatFilter_flag, NucleotideRepeatFilter_threshold):
                    #calculate primer tm
                    primer_tm = float(NN_Tm(seq = primer, compl_seq=complement(primer), primer_conc=primer_conc, \
                                      Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, ion_corr=True))
                    #to determine whether primer tm within a given ranges
                    if primer_tm >= minTm and primer_tm <= maxTm:
                        #calculate form hairpin and homodimer need tm
                        hairpin_tm = HairPinTm(primer, monovalent_cation_eq, primer_conc)
                        homodimer_tm = HomodimerTm(primer, monovalent_cation_eq, primer_conc)
                        #estimate whether can form hainpin or homodimer under the primer tm
                        if primer_tm-float(hairpin_tm) >= float(self_Tmdiff) and primer_tm-float(homodimer_tm) >= float(self_Tmdiff):
                            #write meet condition primer to primer_dict
                            if primer not in primer_dict:
                                primer_dict[primer] = [num, chrom, start, ploc, len(primer), 1, primer_tm, strand]

                            else:
                                primer_dict[primer][4] += 1

'''Function: chop input sequence within the give range 
    Input: 
          @sequence                             : need to amplicon region requence ATGCGTG...[chr7..55018901...55019565]
          @primer_num                           : the region num from bed file
          @chrom,start,minTm,maxTm,GC_range_min,GC_range_max,CheckATends_flag             : follow the functon PrimerFilter 
                                             eg : chr7,55018901,50,70,40,60,1
          @CheckGCclamp_flag,NucleotideRepeatFilter_flag,NucleotideRepeatFilter_threshold : follow the function PrimerFilter
                                             eg : 1,1,4
          @self_Tmdiff, monovalent_cation_eq : follow the function PrimerFilter
                                             eg : 20, global function()
                                                                         
    Output: 
          @primer_dict                          : follow the PrimerFilter function 
                        {"ATGCGCGTGCGTGCATGCGTC": [2,chr7,55018901,100,24,1,70.5,"+"]}
                                                   |  |      |      |   |    |   |
                                                primer_num  region_start len  tm strand
                                                    chrom         primer_start_from_sequence
'''
def ChopImputSeq(sequence, primer_num, primer_size, chrom, start, minTm, maxTm, GC_range_min, GC_range_max, CheckATends_flag, \
                 CheckGCclamp_flag, NucleotideRepeatFilter_flag, NucleotideRepeatFilter_threshold, self_Tmdiff, monovalent_cation_eq):
    primer_dict = {}
    for i in xrange(len(sequence)-primer_size+1):
        primer = sequence[i:i+primer_size].upper()      #positive primer
        rev_comp_primer = rev_complement(primer)        #negative primer

        #primer filter use PrimerFilter function; filter positive and negative primer
        PrimerFilter(primer_num, primer, chrom, start, i, "+", primer_dict, minTm, maxTm,GC_range_min, GC_range_max, CheckATends_flag, \
                     CheckGCclamp_flag, NucleotideRepeatFilter_flag, NucleotideRepeatFilter_threshold, self_Tmdiff, monovalent_cation_eq)
        PrimerFilter(primer_num, rev_comp_primer, chrom, start, i, "-", primer_dict, minTm, maxTm,GC_range_min, GC_range_max, CheckATends_flag, \
                     CheckGCclamp_flag, NucleotideRepeatFilter_flag, NucleotideRepeatFilter_threshold, self_Tmdiff, monovalent_cation_eq)
    return primer_dict


def main():
    seq_file = open(outdir+"sequence.txt", "r")
    #pdb.set_trace()
    parameters_used = open( outdir +'run_summary.txt', 'a')

    unblast_file = open(outdir+"UOD_featureFilter.txt", "w")      #after feature filter
    UOD_all_fasta = outdir+"UOD_featureFilter.fasta"
    FRprimer = open(UOD_all_fasta, "w")                           #similar with unblast_file

    FPrimer = open(outdir+"UOD_forward_primer.fasta", "w")          #final all forward primers
    RPrimer = open(outdir+"UOD_reverse_primer.fasta", "w")          #final all reverse primers
    final_UOD_primer = open(outdir+"UOD_final_primer.txt", "w")        #UOD final primer result include all primers include forward and reverse primers
    final_UOD_all_primer_info = open(outdir+"UOD_final_all_primer_info.fasta", "w")
    final_UOD_all_primer_fasta = open(outdir+"UOD_final_all_primer_fasta.fasta", "w")

    unblast_file.write("chrom\tstart\tend\toccurrence\tsequence\tlen\tstrand\ttm\n")
    monovalent_cation_eq    =    mM_monovalent(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs)
    primer_dict_all = dict()
    primer_num = 0
    
    #open TRS.py result sequence.txt
    for line in seq_file.readlines():
        if line.startswith('>'):
            primer_num += 1
            llist = line.rstrip().split('_')
            chrom, start = llist[1], llist[2]
        else:
            #chop region and filter primer
            primer_l_min, primer_l_max = int(primer_size_range.split("-")[0]), int(primer_size_range.split("-")[1])
            for primer_l in range(primer_l_min, primer_l_max+1, 1):
                primer_dict = ChopImputSeq(line.rstrip(), primer_num, primer_l, chrom, start, minTm, maxTm, GC_range_min, GC_range_max, CheckATends_flag, \
                        CheckGCclamp_flag, NucleotideRepeatFilter_flag, NucleotideRepeatFilter_threshold, self_Tmdiff, monovalent_cation_eq)
                primer_dict_all.update(primer_dict)
    no_primers_designed = len(primer_dict_all)
    sorted_primer_dict = sorted(primer_dict_all.iteritems(), key=itemgetter(1), reverse=False)

    #write meet all condition primers to "UOD_featureFilter.txt" and "UOD_forward_primer.fasta" and "UOD_reverse_primer.fasta"
    for primer_info in sorted_primer_dict:
        primer = primer_info[0]
        primer_region_num = primer_info[1][0]
        chrom = primer_info[1][1]
        start = int(primer_info[1][2]) + primer_info[1][3]
        end = start + primer_info[1][4]
        occu = primer_info[1][5]
        length = primer_info[1][4]
        tm = primer_info[1][6]
        strand = primer_info[1][7]
        unblast_file.write("%d\t%s\t%d\t%d\t%d\t%s\t%d\t%s\t%f\n" %(primer_region_num, chrom, start, end, occu, primer, length, strand, tm))
        FRprimer.write(">TA_" + str(primer_region_num) + "_" + chrom + "_" + str(start) + "_" + str(length) + "_" + str(tm) + "_" + str(strand) + "\n")
        FRprimer.write(primer+"\n")
      
    unblast_file.close()
    FRprimer.close()

    #judge whether the genome file is exist
    if not os.path.exists(genome_fasta):
         sys.stderr.write('\r[*] Please give the hg19.fasta file[hg19.fasta]\n')
         exit(-1)

    os.chdir(refDB_path)
    #create database of balstn
    file_set = {"zmv2all.nin", "zmv2all.nhr", "zmv2all.nsq", "zmv2all.nsi", "zmv2all.nsd", "zmv2all.nog"}
    if set(glob.glob("zmv2all.*")) < file_set:
        f0 = open(os.devnull, 'w')
        sp.call(["makeblastdb","-in","%s" %genome_fasta, "-dbtype","nucl","-parse_seqids","-out","%szmv2all" %refDB_path], stdout=f0,stderr=f0)

    word_size = int(primer_size_range.split("-")[0])      #exact blast word size
    fasta_input_file = UOD_all_fasta
    #make exact balst for the forward primer sequence
    p1 = sp.Popen(["blastn","-task","blastn","-db","%szmv2all" %refDB_path,"-query","%s" %fasta_input_file,"-evalue","%s" %em_e_value,"-word_size","%s"  \
                %word_size,"-gapopen","%s" %em_gapopen,"-gapextend","%s" %em_gapextend,"-reward","%s" %em_reward,"-penalty","%s" %em_penalty,"-dust","no", \
                "-perc_identity","%s" %em_perc_identity,"-max_target_seqs","%s" %em_max_target_seqs,"-max_hsps","%s" %em_max_hsps, \
                "-outfmt","10 qseq qlen qseqid sacc sstart send sstrand", "-num_threads","%s" %em_num_threads],stdout=sp.PIPE)
    exact_match_output, error = p1.communicate()

    #process the blast result
    exact_match_set    =    set() 
    for exact_match_output_line in exact_match_output.split('\n')[:-1]:
        print("%s" %(exact_match_output_line))
        exact_match_output_line = exact_match_output_line.strip(' ').split(',')
        Primer            =    exact_match_output_line[0]
        qseqid            =    exact_match_output_line[2].split('_')
        qseq_chr        =    qseqid[2]
        qseq_start    =    int(qseqid[3])
        qseq_strand    =    qseqid[6]
        qseq_stop    =    int(qseq_start) + int(qseqid[4])
        targetseq_chr    =    exact_match_output_line[3]
        targetseq_start    =    int(exact_match_output_line[4])
        targetseq_stop    =    int(exact_match_output_line[5])
        alignment_length    =    len(Primer)
        query_length    =    int(exact_match_output_line[1])
        if alignment_length == query_length:
            if qseq_chr != targetseq_chr:       #off-target primer sequence//chrom dont same
                exact_match_set.add(Primer)
            if qseq_chr == targetseq_chr:
                if qseq_strand == "+":
                    if (qseq_start+1) != targetseq_start and (qseq_stop-1) != targetseq_stop:   #off-target primer sequence//chrom same but position dont same
                        exact_match_set.add(Primer)
                if qseq_strand == "-":
                    if (qseq_start+1) != targetseq_stop and (qseq_stop-1) != targetseq_start:   #follow the same
                        exact_match_set.add(Primer)

    ### Remove from original dictionary, those primers with exact matches elsewhere in the genome
    for primer_exact_match in exact_match_set:
        if primer_exact_match in primer_dict_all:
            primer_dict_all.pop(primer_exact_match, None)       #pop the primer have the off-target
    no_primers_no_exact_match   =   len(primer_dict_all)    

    #write the  primer sequence after exact balst
    final_UOD_primer.write("chrom\tstart\tend\tseq\ttm\tstrand\n")
    for primer, value in primer_dict_all.items():
        primer_region_num = value[0]
        chrom = value[1]
        primer_start_pos = int(value[2]) + value[3]
        primer_end_pos = int(value[2]) + value[3] + value[4]
        tm = value[6]
        strand = value[7]
        if strand == "+":
            FPrimer.write(">TA_" + str(primer_region_num) + "_" + chrom + "_" + str(primer_start_pos) + "_" + str(len(primer)) + "_" + str(tm)+ "_" + strand + "_" + primer + "\n")
            FPrimer.write(primer + "\n")
        if strand == "-":
            RPrimer.write(">TA_" + str(primer_region_num) + "_" + chrom + "_" + str(primer_start_pos) + "_" + str(len(primer)) +  "_" + str(tm)+ "_" + strand +"_" + primer + "\n")
            RPrimer.write(primer + "\n")
        final_UOD_primer.write(str(primer_region_num) + "\t" + chrom  + "\t" + str(primer_start_pos) + "\t" + str(primer_end_pos) + "\t" + primer + "\t" + str(tm) + "\t" + strand + "\n")
        final_UOD_all_primer_info.write(">TA_" + str(primer_region_num) + "_" + chrom + "_" + str(primer_start_pos) + "_" + str(len(primer)) + "_" + str(tm) + "\n")
        final_UOD_all_primer_info.write(primer + "\n")
        final_UOD_all_primer_fasta.write(">TA_" + str(primer_region_num) + "_" + chrom + "_" + str(primer_start_pos) + "_" + str(len(primer)) + "_" + str(tm)+ "_" + strand + "_" + primer + "\n")
        final_UOD_all_primer_fasta.write(primer + "\n")

    ############################################################
    #Time to run the code: end timer
    ############################################################
    t1 = time.time()
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))
    parameters_used.write   (   "no. of primers designed based on filter criteria : " + str(no_primers_designed)+'\n'
                                "no. of primers without exact match               : " + str(no_primers_no_exact_match)+'\n'
                                "### UOD run duration : " + str(total) + " seconds"+'\n'
                                "##########################################################"+"\n"+
                                "\n"+"\n")
    parameters_used.close()


if __name__ == '__main__':
    main()
