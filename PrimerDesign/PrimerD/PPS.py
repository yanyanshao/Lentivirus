#!usr/bin/python

import os
import csv
import pdb
import time
import primer3
import networkx as nx
from parameters import *
import subprocess as sp
import sequence_processing_functions as spf
from Santalucia_NN_Tm import mM_monovalent

t0 = time.time()

monovalent_cation_eq    =    mM_monovalent(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs)
monovalent_salts = Na + K + (Tris/2)

def seq_extraction_loci(locus, id_r, p_start_pos, p_stop_pos):
    with open(locus) as input_fasta:
        line = input_fasta.read()
        lines = line.split("\n")
        id_reg = lines[(int(id_r)-1)*2]
        sequence = lines[(int(id_r)-1)*2+1]
        reg_start, reg_end = int(id_reg.split('_')[2]), int(id_reg.split('_')[3])
        p_start_pos = max(p_start_pos, reg_start)
        p_stop_pos = min(p_stop_pos, reg_end)
        sequence_reg = sequence[(int(p_start_pos) - int(id_reg.split('_')[2])):(int(p_stop_pos) - int(id_reg.split('_')[2]))]

    return sequence_reg



def heterodimer_dg(seq1, seq2, mv_cation=0,primer_conc=0):
        dg =  (primer3.calcHeterodimer(seq1, seq2,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=60, max_loop=30)).tm
        return float(("{0:.2f}".format(round(dg,2))))

### function to calculate max_cover_sets from all dijkstra's shortest paths
def dijkstra_max_cover_amplicon_sets(dijkstra_shortest_paths):
    max_cover_sets  = []
    if len(dijkstra_shortest_paths) == 1:
        max_cover_sets.append(dijkstra_shortest_paths)
    else:
        for shortpath in dijkstra_shortest_paths:
            shortpath_min, shortpath_max = int(shortpath[0][2]), int(shortpath[-1][-1])
            if not len(max_cover_sets):
                max_cover_sets.append(shortpath)
            else:
                flagGoodToInsert = 1
                for selSet in max_cover_sets:
                    selSetStart, selSetEnd = int(selSet[0][2]), int(selSet[-1][-1])
                    # ignore amplicon sets if the corresponding range has already been covered
                    if shortpath_min > selSetStart and shortpath_max < selSetEnd:
                        flagGoodToInsert = 0
                        break
                    # ignore amplicon sets if the corresponding range has already been covered; same start-stop,more # primers
                    if shortpath_min == selSetStart and shortpath_max == selSetEnd and (len(shortpath) >= len(selSet)): 
                        flagGoodToInsert = 0
                        break
                    # ignore amplicon sets if the corresponding range has already been covered; same start
                    if shortpath_min==selSetStart and shortpath_max< selSetEnd:
                        flagGoodToInsert = 0
                        break
                    # ignore amplicon sets if the corresponding range has already been covered; same stop
                    if shortpath_min > selSetStart and shortpath_max==selSetEnd:
                        flagGoodToInsert = 0
                        break
                if flagGoodToInsert:
                    pairsToBeRemoved = []
                    for selSet in max_cover_sets:
                        selSetStart, selSetEnd = int(selSet[0][2]), int(selSet[-1][-1])
                        if shortpath_min <= selSetStart and shortpath_max >= selSetEnd:
                            pairsToBeRemoved.append(selSet)
                    if len(pairsToBeRemoved):
                        for removePair in pairsToBeRemoved:
                            max_cover_sets.remove(removePair)
                    max_cover_sets.append(shortpath)
    return sorted(max_cover_sets)

### function to retrieve amplicons from each multiplex set:
def multiplx_output_create(list_amplicons, output_file1, primer_info_all, amplicon_info_all):
    for amplicon_coord_info in list_amplicons:
        primer_seq    = primer_info_all[amplicon_coord_info]['primer_sequence'].split("_")[-1]
        primer_name   = primer_info_all[amplicon_coord_info]['primer_name']
        primer_reg, primer_chr, primer_start, primer_end = primer_name.split("_")[1], primer_name.split("_")[2], \
                                            primer_name.split("_")[3], int(primer_name.split("_")[3]) + int(primer_name.split("_")[4])
        if primer_name.split("_")[-1] == "R":
            strand = "-"
        if primer_name.split("_")[-1] == "F":
            strand = "+"
        amplicon_name   = primer_name
        amplicon_seq    = amplicon_info_all[amplicon_coord_info]['amplicon_seq']
        output_file1.write(primer_reg + '\t' + primer_chr + '\t'+ str(primer_start) + '\t' + str(primer_end) + '\t' + strand + "\t" + primer_seq +'\n')
    output_file1.close()

def PickPrimerPair(inputfile):
    G = nx.Graph()

    primerList = []
    d_coords_plus = {}
    d_coords_minus = {}
    primer_pair_counter = 0

    f= open(inputfile)
    locus = outdir + "sequence.txt"
    output_primer_pairs = open(outdir+"PPS_primer_pairs_info.txt", "w")
    output_primer_pairs.write('Primer_pair#'+'\t'+'Primer_name'+'\t'+'Primer_seq'+'\t'+ 'Strand'+'\t'+ '5prime_pos' +'\t'+ '3prime_pos' +'\t'+ 'Tm' +'\t'+ \
    'Max_misprimeTm'+'\t'+ 'GC' +'\t'+ 'Primer_dimer_dG' +'\t'+ 'Amplicon_size' +'\t'+ 'Amplicon_GC' +'\t'+ 'Gaps'+'\t'+ 'Polymorphisms'+'\t'+ 'Amplicon_seq'+'\n')

    order_primer_pairs = open(outdir + "_" +"primer_pairs_order.txt", "w")
    order_primer_pairs .write('Primer_name'+'\t'+'Primer_seq'+'\t'+ 'Tm' +'\n')
    

    csv_f = csv.DictReader(f, delimiter=',')

    for row in csv_f:
        id_p, primer, start, end, strand, Tm, max_misprimeTm, GC , Continuous_GC = row['id'], row['Primer'], int(row['start']), \
                int(row['stop']), row['Strand'], float(row['Primer_Tm']), row['Max_misprime_Tm'], float(row['Primer_GC']), int(row['Continuous_GC'])
        if max_misprimeTm == "3prime_mismatch":
            max_misprimeTm = 0
        else:
            max_misprimeTm = float(max_misprimeTm)
        if Tm - max_misprimeTm >= pair_misprimeTm_diff:
            primerList.append(id_p)
            G.add_node(id_p, p_chrom=id_p.split('_')[2], p_start_pos=int(id_p.split('_')[3]), p_stop_pos=int(id_p.split('_')[3])+int(id_p.split('_')[4]), \
                    strand=strand, Tm=Tm, max_misprimeTm=max_misprimeTm, GC=GC, Continuous_GC=Continuous_GC)

    unique_primers = {}
    primer_pair_coords = []
    amplicon_len_list = []
    amplicon_coords_list = []
    mplex_ampl_coords_list = []

    f8 = open(outdir +'bed_separate_tracks_selected_oligos.bed', 'w')
    #f8.write('browser position chr '+str(chr_no)+':'+ str(start_pos)+'-'+str(stop_pos)+'\n')
    f8.write('track name="Primers" description="Primers on separate tracks" visibility=2 colorByStrand="255,0,0 0,0,255"' + '\n')

    Gfor_subsetting=nx.Graph()

    for (nodeId, data) in G.nodes(data=True):
        '''data    ------ {'p_chrom': 'chr17', 'p_start_pos': 7565168, 'Tm': 65.84, 'GC': 52.17, 'p_stop_pos': 7565191, 'max_misprimeTm': 54.02, 'strand': '+'}
           nodeId  ------'TA_1_chr17_7565168_23_65.84_+_TCCCTGGTTAAGAGATCCTCCTG'
        '''
        if data['strand']   == "+":
            forward_chrom   = data['p_chrom']
            forward_start   = data['p_start_pos']
            forward_stop    = data['p_stop_pos']
            amplicon_start  = int(forward_stop)
            amplicon_end    =  int(forward_start) + 120
            forward_Tm      =   int(data['Tm'])
            forward_GC      =   int(data['GC'])
            forward_MaxMispTm      =   int(data['max_misprimeTm'])
          
            #query_nodeId    =   nodeId
            id_f = nodeId.split('_')[1]
            f_primer_length   =   len(nodeId.split('_')[-1])
            primer_name_f     =   "TA_" + str(id_f) + "_" + forward_chrom + "_" + str(data['p_start_pos']) + "_" + str(f_primer_length) + "_F"
        
            if forward_Tm - forward_MaxMispTm >= pair_misprimeTm_diff:
                amplicon_seq    =    seq_extraction_loci(locus, id_f, amplicon_start, amplicon_end)
                if amplicon_seq == '':
                    continue
                if amplicon_gap_filter ==   1:
                    if 'N'*100 in amplicon_seq:
                        continue

                    amplicon_gc    =    spf.gc_content(amplicon_seq.upper())

                    if 'N'*100 in amplicon_seq:
                        gaps    =   'Yes'
                    else:
                        gaps    =   'No'
                    if 'n' in amplicon_seq:
                        indel    =   'Yes'
                    else:
                        indel    =   'No'

                amplicon_coords    = (id_f, forward_chrom, data['p_stop_pos'], amplicon_end)
                amplicon_coords_list.append(amplicon_coords)

                if data['Tm'] >= multiplex_Tm_min and data['Tm'] <= multiplex_Tm_max and data['Continuous_GC'] <= Continus_GC_num:
                    primer_info   =   {'primer_chr':forward_chrom, 'primer_name':primer_name_f, 'primer_sequence':nodeId, \
                                       'strand':data['strand'], 'p_start_pos':data['p_start_pos'], 'p_stop_pos':data['p_stop_pos'], \
                                       'Tm':data['Tm'], 'max_misprimeTm':data['max_misprimeTm'], 'GC':data['GC']}
                    amplicon_info   =   {'amplicon_gc':amplicon_gc,'gaps':gaps,'indel':indel, 'amplicon_seq':amplicon_seq}

                    Gfor_subsetting.add_node(amplicon_coords, id_f = id_f,primer_info = primer_info, amplicon_info = amplicon_info)
                    mplex_ampl_coords       = (id_f, forward_chrom, data['p_stop_pos'], amplicon_end)
                    mplex_ampl_coords_list.append(mplex_ampl_coords)
            
                    primer_pair_counter += 1
                    output_primer_pairs.write(str(primer_pair_counter) + '\t' + str(primer_name_f)+'\t'+str(nodeId.split("_")[-1])  +'\t'+ \
                                          str(data['strand']) +'\t'+ str(data['p_start_pos']) +'\t'+ str(data['p_stop_pos']) +'\t'+ str(data['Tm'])+'\t'+ \
                                          str(data['max_misprimeTm'])+'\t'+ str(data['GC'])  +'\t'+ \
                                          str(amplicon_gc)+'\t'+ str(gaps)+'\t'+ str(indel)+'\t'+ str(amplicon_seq)+'\n')

        if data['strand']   == "-":
            reverse_chrom   = data['p_chrom']
            reverse_start   = data['p_start_pos']
            reverse_stop    = data['p_stop_pos']
            amplicon_start  = int(reverse_start) - 120
            amplicon_end    = int(reverse_start)
            id_r = nodeId.split('_')[1]

            reverse_Tm      =   int(data['Tm'])
            reverse_GC      =   int(data['GC'])
            reverse_MaxMispTm      =   int(data['max_misprimeTm'])
            r_primer_length   =   len(nodeId.split('_')[-1])
            primer_name_r     =   "TA_" + str(id_r) + "_" +reverse_chrom + "_" + str(data['p_start_pos']) + "_" + str(r_primer_length) + "_R"

            if reverse_Tm - reverse_MaxMispTm >= pair_misprimeTm_diff:
                amplicon_seq    =    seq_extraction_loci(locus, id_r, amplicon_start, amplicon_end)
                if amplicon_seq == '':
                    continue
                if amplicon_gap_filter ==   1:
                    if 'N'*100 in amplicon_seq:
                        continue

                    amplicon_gc    =    spf.gc_content(amplicon_seq.upper())

                    if 'N'*100 in amplicon_seq:
                        gaps    =   'Yes'
                    else:
                        gaps    =   'No'
                    if 'n' in amplicon_seq:
                        indel    =   'Yes'
                    else:
                        indel    =   'No'            

                amplicon_coords    = (id_r, reverse_chrom, amplicon_start, amplicon_end)
                amplicon_coords_list.append(amplicon_coords)

                if data['Tm'] >= multiplex_Tm_min and data['Tm'] <= multiplex_Tm_max and data['Continuous_GC'] <= Continus_GC_num:
                    primer_info   =   {'primer_chr':reverse_chrom, 'primer_name':primer_name_r, 'primer_sequence':nodeId, \
                                       'strand':data['strand'], 'p_start_pos':data['p_start_pos'], 'p_stop_pos':data['p_stop_pos'], \
                                       'Tm':data['Tm'], 'max_misprimeTm':data['max_misprimeTm'], 'GC':data['GC']}
                    amplicon_info   =   {'amplicon_gc':amplicon_gc,'gaps':gaps,'indel':indel, 'amplicon_seq':amplicon_seq}

                    Gfor_subsetting.add_node(amplicon_coords, id_r = id_r,primer_info = primer_info, amplicon_info = amplicon_info)
                    mplex_ampl_coords       = (id_r, reverse_chrom, amplicon_start, data['p_start_pos'])
                    mplex_ampl_coords_list.append(mplex_ampl_coords)

                    primer_pair_counter += 1
                    output_primer_pairs.write(str(primer_pair_counter) + '\t' + str(primer_name_r)+'\t'+str(nodeId.split("_")[-1])  +'\t'+ \
                                          str(data['strand']) +'\t'+ str(data['p_start_pos']) +'\t'+ str(data['p_stop_pos']) +'\t'+ str(data['Tm'])+'\t'+ \
                                          str(data['max_misprimeTm'])+'\t'+ str(data['GC'])  +'\t'+ \
                                          str(amplicon_gc)+'\t'+ str(gaps)+'\t'+ str(indel)+'\t'+ str(amplicon_seq)+'\n')


    
        #primer_length   =   len(nodeId.split('_')[-1])

    f8.close()
    output_primer_pairs.close()
    order_primer_pairs.close()
    no_unique_primers_picked    =   len(unique_primers)

    return {'mplex_ampl_coords_list': mplex_ampl_coords_list, 'Gfor_subsetting': Gfor_subsetting, 'amplicon_coords_list':amplicon_coords_list,'primer_pair_counter': primer_pair_counter}


 

### function to retrieve multiplex compatible groups :
def multiplex_group_output(input_file, set_no, primer_pair_counter, multiplx_pooled_out):
    with open(input_file) as f:
        content = f.readlines()
        for line in content:
            if line.startswith('Group'):
                line = line.split()
                group = line[1]
                for amplicon in line[2:]:
                    chrom, f_primer, r_primer =  amplicon.split('_')[1],amplicon.split('_')[2], amplicon.split('_')[6]
                    f_primer_len,r_primer_len = amplicon.split('_')[3], amplicon.split('_')[7]
                    forward_primer_start = int(f_primer.split('_')[2])
                    for key, value in r_primer_info_all.iteritems():
                        if value['primer_name'] == r_primer and int(key[0]) == forward_primer_start :
                            r_primer_info = value
                            selected_coords = key
                    primer_pair_counter += 1
                    f_primer_name, f_primer_seq, fStrand, fStart_pos, fStop_pos, fTm, fMax_misprimeTm, fGC, Primer_dimer_dG, \
                            Amplicon_size, Amplicon_gc, Gaps, Polymorphisms, Amplicon_seq \
                             = f_primer_info_all[selected_coords]['primer_name'],f_primer_info_all[selected_coords]['primer_sequence'],\
                             f_primer_info_all[selected_coords]['strand'],f_primer_info_all[selected_coords]['p_start_pos'],\
                             f_primer_info_all[selected_coords]['p_stop_pos'],f_primer_info_all[selected_coords]['Tm'],\
                             f_primer_info_all[selected_coords]['max_misprimeTm'],f_primer_info_all[selected_coords]['GC'],amplicon_info_all[selected_coords]['interaction_dg'],\
                             amplicon_info_all[selected_coords]['amplicon_len'],amplicon_info_all[selected_coords]['amplicon_gc'],amplicon_info_all[selected_coords]['gaps'],\
                             amplicon_info_all[selected_coords]['indel'],amplicon_info_all[selected_coords]['amplicon_seq']

                    multiplx_pooled_out.write(str(set_no) +'\t'+ str(group) +'\t'+ str(primer_pair_counter) +'\t'+ f_primer_name+'\t'+ f_primer_seq +'\t'+ \
                            fStrand +'\t'+ str(fStart_pos) +'\t'+ str(fStop_pos) +'\t'+ str(fTm) +'\t'+ str(fMax_misprimeTm) +'\t'+ str(fGC) +'\t'+ \
                            str(Primer_dimer_dG) +'\t'+ str(Amplicon_size) +'\t'+ str(Amplicon_gc) +'\t'+ str(Gaps) +'\t'+ Polymorphisms +'\t'+ Amplicon_seq +'\n')

                    r_primer_name, r_primer_seq, rStrand, rStart_pos, rStop_pos, rTm, rMax_misprimeTm, rGC \
                                    = r_primer_info['primer_name'],r_primer_info['primer_sequence'],r_primer_info['strand'],\
                                    r_primer_info['p_start_pos'],r_primer_info['p_stop_pos'],r_primer_info['Tm'],\
                                    r_primer_info['max_misprimeTm'],r_primer_info['GC']
                    multiplx_pooled_out.write(str(set_no) +'\t'+ str(group) +'\t'+ str(primer_pair_counter) +'\t'+ r_primer_name+'\t'+ r_primer_seq+'\t'+ rStrand +'\t'+ \
                            str(rStop_pos) +'\t'+ str(rStart_pos) +'\t'+ str(rTm) +'\t'+ str(rMax_misprimeTm) +'\t'+ str(rGC) +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\t'+ '-' +'\n')
    f.close()
    return int(primer_pair_counter)

def FindOverLap(filter_primer_list, primer):
    '''filter_primer_list ----[(chr,start,end,strand,sequence)]
       primer----(chr,start,end,strand,sequence)
    '''
    s1, e1, s2, e2 = int(filter_primer_list[-1][1]), int(filter_primer_list[-1][2]), int(primer[1]), int(primer[2])
    
    if s1 <= s2 and s2 <= e1 :
        return True

    return False
    
    

def main():
    #pdb.set_trace()
    input_file = outdir + "PSE_final_result.csv"
    final_primer_file_1 = outdir + "PPS_multiplx_out_pool_1.txt"        #final result file
    final_primer_file_2 = outdir + "PPS_multiplx_out_pool_2.txt"
    pool1  = open(final_primer_file_1, "w")
    pool2  = open(final_primer_file_2, "w")


    ### Parameter out put file
    parameters_used     = open(outdir +'run_summary.txt', 'a')

    ### print parameters used
    parameters_used.write(  "##########################################################"+"\n"+
                    "### Summary of PPS parameters"+"\n"+
                    "### Start date: " + str(time.strftime("%m/%d/%Y"))+";"+ "\t"+ "Start time: "+ str(time.strftime("%H:%M:%S"))+"\n"+
                    "##########################################################"+"\n"+
                    "Max_primer_pair_Tm_diff    =   "+str(pair_Tm_diff)+"\n"+
                    "Min_primer_misprimeTm_diff =   "+str(pair_misprimeTm_diff)+"\n"+
                    "Min_Amplicon_size          =   "+str(amplicon_size_min)+"\n"+
                    "Max_Amplicon_size          =   "+str(amplicon_size_max)+"\n"+
                    "amplicon_gap_filter        =   "+str(amplicon_gap_filter)+ "(1:Yes; 0:No)" +"\n"+
                    "min_misprime_dg            =   "+str(min_misprime_dg)+"\n"+
                    "multiplex_Tm_min           =   "+str(multiplex_Tm_min)+"\n"+
                    "multiplex_Tm_max           =   "+str(multiplex_Tm_max)+"\n"+
                    "##########################################################"+"\n"
                )

    primer_picking_output = PickPrimerPair(input_file)
    mplex_ampl_coords_list, Gfor_subsetting,  amplicon_coords_list, primer_pair_counter = primer_picking_output['mplex_ampl_coords_list'], \
                                                primer_picking_output['Gfor_subsetting'], primer_picking_output['amplicon_coords_list'],  \
                                                primer_picking_output['primer_pair_counter']
   
    sorted_amplicon_coords_list =   sorted(list(set(mplex_ampl_coords_list)))

    for i in xrange(len(sorted_amplicon_coords_list) -1):
        query_node  = sorted_amplicon_coords_list[i]
        query_reg_id = query_node[0]
        query_chrom, query_amplicon_start, query_amplicon_stop = query_node[1], int(query_node[2]), int(query_node[3])

        for j in xrange(len(sorted_amplicon_coords_list)-(i+1)):
            target_node  = sorted_amplicon_coords_list[i+j+1]
            target_reg_id = target_node[0] 
            target_chrom, target_amplicon_start, target_amplicon_stop = target_node[1], int(target_node[2]), int

            if (query_reg_id == target_reg_id) and (query_amplicon_stop >= target_amplicon_start) and \
                    (target_amplicon_stop > query_amplicon_stop) and (target_amplicon_start > query_amplicon_start):
                Gfor_subsetting.add_edge(query_node,target_node)
            if query_reg_id != target_reg_id:
                break
    ### pick all sub networks
    graphs = list(nx.connected_component_subgraphs(Gfor_subsetting))
    bestshortest_paths_list  =   []
    all_short_paths = []

    primer_info_all   =  nx.get_node_attributes(Gfor_subsetting,'primer_info')
    amplicon_info_all   =  nx.get_node_attributes(Gfor_subsetting,'amplicon_info')

    all_nodes_connectedcomponents = []
    graph_no    =   0
    subnetwork_stats = {}
    counter = 1
    subnetwork_stats_included = []

    candidate_dict_pos = {}
    for amplicon, primer in primer_info_all.items():
        '''
        ('1', 'lentivirus', 414, 534): {'max_misprimeTm': 49.78, 'primer_name': 'TA_1_lentivirus_534_29_R', \
                'primer_chr': 'lentivirus', 'GC': 41.38, 'p_stop_pos': 563, 'primer_sequence': 'TA_1_lentivirus_534_29_65.75_-_TGTCTACAGCCTTCTGATGTTTCTAACAG', \
                'Tm': 65.75, 'p_start_pos': 534, 'strand': '-'}
        '''
        p_start, p_end, p_seq, p_strand, p_gc, p_tm, p_len = primer['p_start_pos'], primer['p_stop_pos'], primer['primer_sequence'].split("_")[-1], \
                               primer['primer_sequence'].split("_")[-2], primer['GC'], primer['Tm'], primer['primer_sequence'].split("_")[-4]
        coord = str(p_start) + "-" + str(p_end)

        if p_start <= 2000:
            if p_strand == "-":
                if coord not in candidate_dict_pos:
                    candidate_dict_pos[coord] = [p_seq, p_len, p_strand, p_tm, p_gc]
                else:
                    print("Error may Occur!")
        if p_start >= 3000:
            if p_strand == "+":
                if coord not in candidate_dict_pos:
                    candidate_dict_pos[coord] = [p_seq, p_len, p_strand, p_tm, p_gc]
                else:
                    print("Error may Occur!")


    candidate_list_pos = sorted(candidate_dict_pos.keys(),key=lambda student : int(student.split("-")[0]))
    first_n0 = candidate_list_pos[0]
    pair_20_bp = [first_n0]
    flag_20, keep_p = 1, []
    final_list_pos = []
    for i in xrange(1, len(candidate_list_pos)):
        if int(candidate_list_pos[i].split("-")[0]) - int(first_n0.split("-")[0]) < 20:
            pair_20_bp.append(candidate_list_pos[i])

        if int(candidate_list_pos[i].split("-")[0]) - int(first_n0.split("-")[0]) >= 20 or i == len(candidate_list_pos) - 1:
            if flag_20:
                final_list_pos.append(pair_20_bp[len(pair_20_bp)/2])
                keep_p  = pair_20_bp[len(pair_20_bp)/2]
                flag_20 -= 1
        if keep_p:
            if int(candidate_list_pos[i].split("-")[0]) - int(keep_p.split("-")[0]) > 50:
                first_n0 = candidate_list_pos[i]
                pair_20_bp = [first_n0]
                flag_20 = 1
                keep_p = ''

    pool1.write("lentivirus" + "\t" + "start" + "\t" + "end" + "\t" + "sequence" + "\t" + "length" + "\t" + "strand" + "\t" + "tm" + "\t" + "gc" + "\n")
    pool2.write("lentivirus" + "\t" + "start" + "\t" + "end" + "\t" + "sequence" + "\t" + "length" + "\t" + "strand" + "\t" + "tm" + "\t" + "gc" + "\n")
    n = 0
    for i in final_list_pos:
        n += 1
        info = candidate_dict_pos[i]
        if n%2:
            pool1.write("lentivirus" + "\t" + i.split("-")[0] + "\t" + i.split("-")[1] + "\t" + info[0] + "\t" + str(info[1]) + "\t" + info[2] + "\t" + str(info[3]) + "\t" + str(info[4]) + "\n")
        else:
            pool2.write("lentivirus" + "\t" + i.split("-")[0] + "\t" + i.split("-")[1] + "\t" + info[0] + "\t" + str(info[1]) + "\t" + info[2] + "\t" + str(info[3]) + "\t" + str(info[4]) + "\n")
    

    '''
    for graph in graphs:
        amplicon_coords = sorted(nx.nodes(graph))
        #PosExtract(amplicon_coords)
        amplicon_coords_sorted_by_second = sorted(amplicon_coords, key=lambda tup: tup[3])
        all_nodes_connectedcomponents.append(amplicon_coords)
        start_node  =  amplicon_coords[0]
        stop_node   = amplicon_coords_sorted_by_second[-1]

         ### record the start stop pos for each subnetwork
        if len(amplicon_coords) > 1:
            graph_no    += 1
            subnetwork_stats['subnetwork# '+ str(graph_no)]  = ['chrom = ' + start_node[1], 'start = ' + str(start_node[2]), 'stop = ' + str(stop_node[3]), 'length = '+str(int(stop_node[3]) - int(start_node[2]) + 1)+ ' bp', '# primer pairs = ' + str(len(amplicon_coords))]
            subnetwork_stats_included.append((start_node[0], start_node[1], start_node[2], stop_node[3]))

        max_coverage_len = 0
        DG=nx.DiGraph()
        DG.add_nodes_from(amplicon_coords)

        # first iteration to get max edge weight (penalized cumulative coverage)
        for i in xrange(len(amplicon_coords)-1):
                reg_n0, chrom_n0, start_pos_n0,stop_pos_n0 = amplicon_coords[i][0], amplicon_coords[i][1], amplicon_coords[i][2], amplicon_coords[i][3]
                for j in xrange(len(amplicon_coords[i+1:])):
                    reg_n1, chrom_n1, start_pos_n1,stop_pos_n1 = amplicon_coords[i+j+1][0], amplicon_coords[i+j+1][1], amplicon_coords[i+j+1][2], amplicon_coords[i+j+1][3]
                    if (reg_n0 == reg_n1) and (chrom_n0 == chrom_n1) and(stop_pos_n0 >= start_pos_n1) and (stop_pos_n1 > stop_pos_n0) and (start_pos_n1 > start_pos_n0):
                        coverage    = ((stop_pos_n1 - start_pos_n0)+1) - ((stop_pos_n0 - start_pos_n1) +1)
                        if coverage > max_coverage_len:
                            max_coverage_len = coverage
                    if (reg_n0 != reg_n1) or (chrom_n0 != chrom_n1):
                        break

        
        # second iteration with adjusted edge weight (penalized cumulative coverage) and draw edged with the adjusted edge weight                
        for i in xrange(len(amplicon_coords)-1):
                reg_n0, chrom_n0, start_pos_n0,stop_pos_n0 = amplicon_coords[i][0],amplicon_coords[i][1],amplicon_coords[i][2], amplicon_coords[i][3]
                for j in xrange(len(amplicon_coords[i+1:])):
                    reg_n1, chrom_n1, start_pos_n1,stop_pos_n1 = amplicon_coords[i+j+1][0], amplicon_coords[i+j+1][1],amplicon_coords[i+j+1][2], amplicon_coords[i+j+1][3]
                    if (reg_n0 == reg_n1) and (chrom_n0 == chrom_n1) and (stop_pos_n0 >= start_pos_n1) and (stop_pos_n1 > stop_pos_n0 ) and (start_pos_n1 > start_pos_n0):
                        coverage    = max_coverage_len - (((stop_pos_n1 - start_pos_n0)+1) - ((stop_pos_n0 - start_pos_n1) +1))
                        DG.add_weighted_edges_from([(amplicon_coords[i],amplicon_coords[i+j+1],coverage)])
        try:
            best_shortest_path    =    nx.dijkstra_path(DG,start_node, stop_node)
            #bestshortest_paths_list.append(best_shortest_path)
        except:
            pass

        bestshortest_paths_list.append(best_shortest_path)
    bestshortest_paths_list = dijkstra_max_cover_amplicon_sets(bestshortest_paths_list)
    
        
    all_nodes_connectedcomponents   = [item for sublist in all_nodes_connectedcomponents for item in sublist]
    bestshortest_paths_list         = [item for sublist in bestshortest_paths_list for item in sublist]
    
    
    ### separate odd and even elements of bestshortest_paths_list so that overlapping amplicons are not multiplexed together
    bestshortest_paths_list = sorted(bestshortest_paths_list)

    multiplx_output_create(bestshortest_paths_list, output_file1, primer_info_all, amplicon_info_all)

    all_primer_dict, bed_list = {}, []
    output_file1.close()
    all_primer = open(final_primer_file1, "r")
    bed_target = open("/home/yys/Project/Lilab/Primer/PrimerDesign_SB/test/TP53.bed", "r")

    for line in bed_target.readlines():
        lline = line.rstrip().split("\t")
        bed_list.append((lline[0], lline[1], lline[2]))

    for line in all_primer.xreadlines():
        lline = line.rstrip().split("\t")
        reg, chrom, start, end, strand, seq = lline[0], lline[1], lline[2], lline[3], lline[4], lline[5]
        if reg not in all_primer_dict:
            all_primer_dict[reg] = [(chrom, start, end, strand, seq)]
        else:
            all_primer_dict[reg].append((chrom, start, end, strand, seq))

    final_primer_dict = {}
    for reg, primers in all_primer_dict.items():
        primers = sorted(primers)
        rstart, rend = bed_list[int(reg)-1][1], bed_list[int(reg)-1][2]
        final_primer_list = []
        l, i = len(primers), 0
        
        while(i<l):
            if primers[i][3] == "-" and (int(rstart) - int(primers[i][2]) >=0 or  int(primers[i][1]) - int(rstart) <10):
                i += 1
                continue
            if primers[i][3] == "+" and (int(rend) - int(primers[i][2]) <= 10 or int(primers[i][1]) - int(rend) >= 0 ):
                i += 1
                continue

            filter_primer_list = [primers[i]]
            i += 1
            if i < len(primers):
                while(FindOverLap(filter_primer_list, primers[i])):
                    filter_primer_list.append(primers[i])
                    i += 1
                    if i > len(primers) - 1:
                        break

            if len(filter_primer_list) == 1 :
                final_primer_list.extend(filter_primer_list)
            else:
                flag = 0
                if i-len(filter_primer_list) <= 0:
                    r_primer = primers[i]
                    for p in filter_primer_list:
                        if p[3] == r_primer[3]:
                            flag = 1
                            final_primer_list.append(p)
                            break
                    if flag == 0:
                        final_primer_list.append(filter_primer_list[-1])
                    continue
                if i == len(primers):
                    f_primer = primers[i-len(filter_primer_list)-1]
                    for p in filter_primer_list:
                        if p[3] == f_primer[3]:
                            flag = 1
                            final_primer_list.append(p)
                            break
                    if flag == 0:
                        final_primer_list.append(filter_primer_list[0])
                    continue   
                else:
                    f_primer, r_primer =  primers[i-len(filter_primer_list)-1], primers[i]
                    if int(f_primer[2])-int(filter_primer_list[0][1]) > int(filter_primer_list[-1][2])-int(r_primer[1]):
                        for p in filter_primer_list:
                            if p[3] == f_primer[3]:
                                flag = 1
                                final_primer_list.append(p)
                                break
                        if flag == 0:
                            final_primer_list.append(filter_primer_list[0])
                          

                    else:
                        for p in filter_primer_list:
                            if p[3] == r_primer[3]:
                                flag = 1
                                final_primer_list.append(p)
                                break
                        if flag == 0:
                            final_primer_list.append(filter_primer_list[-1])

        final_primer_dict[reg] = final_primer_list

    for r,ps in final_primer_dict.items():
        for p in ps:
            output_file2.write(p[0] + "\t" + p[1] + "\t" + p[2] + "\t" + p[3] + "\t" + p[4] + "\n")



    ############################################################
    # Time to run the code: end timer
    ############################################################
    t1 = time.time()
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))
    parameters_used.write(
                                "subnetwork coverage stats                    : " + str(subnetwork_stats) + '\n'
                                "# minimum tiling primers                     : " + str(primer_pair_counter) + '\n'
                                "### PPS run duration           : " + str(total) + " seconds"+'\n'
                                "##########################################################"+"\n"
                        )
    parameters_used.close()
    '''

if __name__ == '__main__':
    main()
