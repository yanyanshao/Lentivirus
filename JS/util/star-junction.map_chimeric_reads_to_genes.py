# -*- coding:utf-8 -*-
###################################################################
# File Name: star-junction.map_chimeric_reads_to_genes.py
# Author: yys
# mail: shayy0919@163.com
# Created Time: 2019年07月31日 星期三 15时30分07秒
###################################################################
#!/usr/bin/python3

import os
import sys
import pdb
from optparse import OptionParser

CHR= {"chr1":249250621, "chr2":243199373, "chr3":198022430, "chr4":191154276, \
       "chr5":180915260, "chr6":171115067, "chr7":159138663, "chr8":146364022, \
       "chr9":141213431, "chr10":135534747, "chr11":135006516, "chr12":133851895,\
       "chr13":115169878, "chr14":107349540, "chr15":102531392, "chr16":90354753,\
       "chr17":81195210, "chr18":78077248, "chr19":59128983, "chr20":63025520, \
       "chr21":48129895, "chr22":51304566, "chrX":155270560, "chrY":59373566, "chrM":16569}


class Node(object):
    def __init__(self, left, right, p, inter, name, maxx):
        self.key = int(inter.low)
        self.left = left
        self.right = right
        self.p = p
        self.name = name
        self.inter = inter  #interval info
        self.maxx = int(maxx)


class Inter(object):
    def __init__(self, low, high, name):
        self.low = int(low)
        self.high = int(high)
        self.name = name


class Tree(object):
    def __init__(self, root, nil):
        self.root = root
        self.nil = nil

    def insert_tree(self, z):
        y = self.nil
        x = self.root

        while x != self.nil:
            y = x
            if z.key < x.key:
                x = x.left
            else:
                x = x.right

        z.p = y
        if y == self.nil:
            self.root == z
        else:
            if z.key < y.key:
                y.left = z
            else:
                y.right = z

        z.left = self.nil
        z.right = self.nil
        z.maxx = max(z.inter.high,z.left.maxx,z.right.maxx)
        
        while(z.p != self.nil):
            z.p.maxx= max(z.p.maxx, z.maxx)
            z = z.p


    def search_all(self, pos):
        x = self.root
        while x!=self.nil and not (pos<=x.inter.high and pos>=x.inter.low):
            if x.left.maxx > pos:
                x = x.left
            else:
                x = x.right

        return x

    def print_tree(self, z):
        if z is None:
            return

        if z != None:
            self.print_tree(z.left)
            print("%d\t%s\t%s" %(z.inter.low, z.inter.high, z.name))
            self.print_tree(z.right)


def BuildIntervalTree(gtffeaturefile):
    gtfp = open(gtffeaturefile, "r")

    gtfdict = {}
    prechr = 'dummy'

    l = 0
    inter=Inter(0,0,'')
    nil=Node(None,None,None,inter,'',0)
    for line in gtfp.readlines():
        l += 1
        sys.stdout.write('%d\r' % l)
        sys.stdout.flush()
        llist = line.rstrip().split()

        if llist[0] == prechr:
            name = llist[-1]
            inter = Inter(int(llist[3]), int(llist[4]), name)
            z=Node(nil, nil, nil, inter, name, 0)
            T.insert_tree(z)
        else:
            if prechr != 'dummy':
                gtfdict[prechr] = T

            prechr = llist[0]
            inter=Inter(CHR[prechr]/2, CHR[prechr]/2, '')
            root=Node(nil,nil,nil,inter,"NONE",CHR[prechr]/2)
            T=Tree(root,nil)
    gtfdict[prechr] = T

    return gtfdict


def BuildGeneTree(genefeaturefile):
    genefp = open(genefeaturefile, "r")

    genedict = {}

    n = 0
    inter=Inter(0,0,'')
    nil=Node(None,None,None,inter,'',0)
    for line in genefp.readlines():
        lline = line.rstrip().split()
        
        n += 1
        sys.stdout.write('%d\r' % n)
        sys.stdout.flush()
        name = lline[5]
        inter = Inter(int(lline[2]), int(lline[3]), name)
        if lline[1] in genedict:
            z=Node(nil, nil, nil, inter, name, 0)
            genedict[lline[1]].insert_tree(z)
        else:
            name = lline[5]
            inter=Inter(int(lline[2]), int(lline[3]), name)
            root=Node(nil,nil,nil,inter,name,0)
            T=Tree(root,nil)
            genedict[lline[1]] = T

    return genedict


def main():
    usage = "\n\tstar-junction.map_chimeric_reads_to_genes.py [options] <ampliconfile.csv> <Chimeric.out.junction> <Aligned.out.bam>\n\n \
    \t###########################################################################\n \
    \t#\n \
    \t# Required:\n \
    \t#     --JuncM | string   Junction file from star-fusion.Match_length.py script[xxx.Mathch_length]\n \
    \t#     --Gene  | string   Gene region file[ref_annot.gtf.gene_spans]\n \
    \t#     --Exon  | string   Exon region file[ref_annot.gtf.mini.sortu.tmp]\n \
    \t#      gene file and exon file from \"/home/yys/Project/Lilab/Lentivirus/DataBase/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/\"\n \
    \t#\n \
    \t############################################################################"

    parser = OptionParser(usage=usage)
    parser.add_option('--JuncM', action='store', dest='JuncM',
            help="[required] Junction file from star-fusion.Match_length.py script[xxx.Mathch_length]")
    parser.add_option('--Gene',  action='store', dest='Gene',
            help="[required] Gene region file[ref_annot.gtf.gene_spans]")
    parser.add_option('--Exon',  action='store', dest='Exon',
            help="[required] Exon region file[ref_annot.gtf.mini.sortu.tmp]")

    (option, args) = parser.parse_args()

    if (not option.JuncM) or (not option.Gene) or (not option.Exon):
        parser.print_help()
        sys.exit(0)

    #BuildIntervalTree("/home/yys/Project/Lilab/Lentivirus/DataBase/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/ref_annot.gtf.mini.sortu.tmp")

    junc = open(option.JuncM)

    gene = BuildGeneTree(option.Gene)
    exon = BuildIntervalTree(option.Exon)

    failure_pass = open(os.path.join(os.path.split(option.JuncM)[0],"star-fusion.junction_breakpts_to_genes_fail.txt"), "w")
    success_pass = open(os.path.join(os.path.split(option.JuncM)[0],"star-fusion.junction_breakpts_to_genes_pass.txt"), "w")


    inter=Inter(0,0,'')
    nil=Node(None,None,None,inter,'',0)
    for line in junc.readlines():
        lline = line.rstrip().split()

        if lline[0] == "chrLT":
            T_exon = exon[lline[3]]
            T_gene = gene[lline[3]]
            coor_exon = T_exon.search_all(int(lline[4]))
            coor_gene = T_gene.search_all(int(lline[4]))
        else:
            T_exon = exon[lline[0]]
            T_gene = gene[lline[0]]
            coor_exon = T_exon.search_all(int(lline[1]))
            coor_gene = T_gene.search_all(int(lline[1]))

        if not coor_gene.name:
            failure_pass.write("%s" %line)
        if coor_gene.name and not coor_exon.name:
            success_pass.write("%s\t%s\tOnly Gene\n" %(line.rstrip(),coor_gene.name))
        if coor_exon.name:
            success_pass.write("%s\t%s\tExon\n" %(line.rstrip(),coor_exon.name))



if __name__ == '__main__':
    main()
