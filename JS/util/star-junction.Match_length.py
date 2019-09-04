# -*- coding:utf-8 -*-

# Version:1.0

###################################################################
# File Name: star-junction.Match_length.py
# Author: yys
# mail: shayy0919@163.com
# Created Time: 2019年07月30日 星期二 08时33分34秒
###################################################################
#!/usr/bin/python3

import os
import re
import sys
import pysam
from optparse import OptionParser

def MaxPrimerLen(filename):
    '''filename == "/home/yys/Project/Lilab/Lentivirus/Data/SecondPrimData/LG.primer"
       #Chr        Start   End          Seq                Strand
       lentivirus   121    148 TTGAAGCACTCAAGGCAAGCTTTATTG   -
    '''
    fp = open(filename, "r")

    lmax = 0
    for line in fp.readlines():
        if line[0] == "#":
            continue

        lline = line.rstrip().split()
        lmax = len(lline[3]) if (len(lline[3])>lmax) else lmax

    return lmax


def GetGenomeCoordsViaCigar(rst, cigar):
    '''rst : column 11 & 13
       cigar : column 12 & 14
       genome_coords : list;match base position in genome
       query_coords : list;match base position in read
       rst  : 68175467      genome_coords :[[68175467,68175529]]
       cigar: 24S62M        query_coords : [[24,62]]
    '''
    genome_lend, query_lend = int(rst), 0
    genome_coords, query_coords = [], []
    genome_lend -= 1

    for ci in re.findall('(\d+)([A-Zp])', cigar):
        mlen, code = int(ci[0]), ci[1]

        if not re.findall('^[MSDNIHp]$', code):
            sys.stderr.write("Error, cannot parse cigar code [%s] " %(code))
            sys.exit(-1)

        if (code == "M"):
            genome_rend = genome_lend + mlen
            query_rend = query_lend + mlen

            genome_coords.append([genome_lend,genome_rend])
            query_coords.append([query_lend, query_rend])

            genome_lend = genome_rend
            query_lend = query_rend

        if (code == "D" or code == "N" or code == "p"):
            genome_lend += mlen

        if (code == "I" or code == "S" or code == "H" ):
            query_lend += mlen

    return (genome_coords, query_coords)
            

def ComputeAnchorLength(read_coords_aref, left_or_right, orient):
    if (left_or_right == "left"):
        if orient == "+":
            term_segment = read_coords_aref[-1]
        else:
            term_segment = read_coords_aref[0]
    else:
        if orient == "+":
            term_segment = read_coords_aref[0]
        else:
            term_segment = read_coords_aref[-1]

    seg_length = term_segment[1] - term_segment[0]

    return seg_length


def main():
    usage = "\n\tstar-fusion.Match_length.py [options] <ampliconfile.csv> <Chimeric.out.junction> <Aligned.out.bam>\n\n \
    \t###########################################################################\n \
    \t#\n \
    \t# Required:\n \
    \t#     --Amplicon | string   Primer Amplicon file[tab]\n \
    \t#     --Chimeric | string   Chimeric.out.junction file\n \
    \t#     --Alignbam | string   Aligned.out.bam file\n \
    \t############################################################################"

    parser = OptionParser(usage=usage)
    parser.add_option('--Amplicon', action='store', dest='Amp',
            help="[required] input primer amplicon file")
    parser.add_option('--Chimeric',  action='store', dest='Chimeric',
            help="[required] junction of Chimeric.out.junction file")
    parser.add_option('--Alignbam',  action='store', dest='Align',
            help="[required] input bam file of STAR result.[Aligned.out.bam]")

    (option, args) = parser.parse_args()

    if (not option.Amp) or (not option.Chimeric):
        parser.print_help()
        sys.exit(0)

    chimericOut = open(os.path.join(os.path.split(option.Chimeric)[0], "std.Chimeric.out.junction.Mathch_length"), "w")
    alignIn = pysam.AlignmentFile(option.Align, "rb")
    chimericIn = open(option.Chimeric, "r")
    pmax = MaxPrimerLen(option.Amp)

    bamDict = {}
    for read in alignIn.fetch(until_eof = True):
        if not read.is_unmapped:
            if read.qname not in bamDict:
                bamDict[read.qname] = 1
            else:
                continue

    for line in chimericIn.readlines():
        lline = line.rstrip().split()
        if lline[9] in bamDict:
            continue
        A_rst, A_cigar, A_orient = lline[10], lline[11], lline[2]
        B_rst, B_cigar, B_orient = lline[12], lline[13], lline[5]

        A_read_coords_aref = GetGenomeCoordsViaCigar(A_rst, A_cigar)[1]
        B_read_coords_aref = GetGenomeCoordsViaCigar(B_rst, B_cigar)[1]

        A_seg_length = ComputeAnchorLength(A_read_coords_aref, "left", A_orient)
        B_seg_length = ComputeAnchorLength(B_read_coords_aref, "right", B_orient)

        if (A_seg_length > pmax and B_seg_length > pmax and \
                lline[0] != lline[3] and (lline[0] == "chrLT" or lline[3] == "chrLT")):
            chimericOut.write("%s\tleft_anchor_length:%d\tright_anchor_length:%d\n" %(line.rstrip(), A_seg_length, B_seg_length))


if __name__ == '__main__':
    main()

