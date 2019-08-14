################################################################
#********************    TRS.py    ****************************#
###############################################################

refDB_path = "/home/yys/Project/Lilab/Primer/PrimerDesign_lentivirus_SB/Database/Database_Venus/"
genome_fasta = "/home/yys/Project/Lilab/Primer/PrimerDesign_lentivirus_SB/Database/Database_Venus/hg19.fasta"
prog = "/home/yys/Project/Lilab/Primer/PrimerDesign_lentivirus_SB/PrimerD/"
cart_fa = "/home/yys/Project/Lilab/Primer/PrimerDesign_lentivirus_SB/LentiGenome/pLenti-EF1-Venus.fasta"
outdir = "/home/yys/Project/Lilab/Primer/PrimerDesign_lentivirus_SB/Result/Venus/"

flanking_flag           =   0       #whether to extend from bed file region
flanking_size           =   50     #extend size


################################################################
#********************    UOD.py    ****************************#
################################################################
primer_size_range   =   "22-35"     #primer size
minTm = 65                      #min primer tm
maxTm = 73      #max primer tm
GC_range_min = 40   #min GC content
GC_range_max = 60   #max GC content
CheckATends_flag = 0    #whether to check AT end
CheckGCclamp_flag = 1   #whether to check GC clamp
Continus_GC_num = 3
NucleotideRepeatFilter_flag = 1     #whether to check repeat sequence
NucleotideRepeatFilter_threshold = 4 #the candidate with repeat > threshold will be discarded
self_Tmdiff = 20

primer_conc = 100              #primer concentration
Na                  =   0      #mM
K                   =   50     #mM
Tris                =   10     #mM
Mg                  =   1.5    #mM  
dNTPs               =   0.2    #mM


### blastn parameters for exact-match(em) search
em_e_value          =   30000
em_gapopen          =   2
em_gapextend        =   2
em_reward           =   1       #reward
em_penalty          =   -3      #penalty
em_perc_identity    =   100     #similarity 100;exact match
em_max_target_seqs  =   2       #find max target sequence num
em_max_hsps         =   2
em_num_threads      =   20       #threads num


################################################################
#********************    PSE.py    ****************************#
################################################################
### blastn parameters for mis-prime(mp) search
mp_e_value          =   30000
mp_gapopen          =   2
mp_gapextend        =   2
mp_reward           =   1   #reward score
mp_penalty          =   -1  #penalty score
mp_num_threads      =   30
mp_perc_identity    =   70      #similarity 70
mp_max_target_seqs  =   3       #find max target sequence num
mp_max_hsps         =   20


################################################################
#********************    PPS.py    ****************************#
###############################################################
amplicon_size_min   =   90     #min amplicon size
amplicon_size_max   =   150     #max amplicon size
amplicon_gap_filter =   1
pair_misprimeTm_diff=   10      #mismtch primer(off-target primer) diff with target-primer
pair_Tm_diff        =   10      # pair primer tm diff
pair_GC_diff        =   20      # pair primer GC diff
min_misprime_dg     =   -2000

multiplex_Tm_min    =   55
multiplex_Tm_max    =   70
