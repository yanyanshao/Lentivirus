# -*- coding:utf-8 -*-
###################################################################
# File Name: star-junction.fiter.py
# Author: yys
# mail: shayy0919@163.com
# Created Time: 2019年08月05日 星期一 07时32分36秒
###################################################################
#!/usr/bin/python3

import os
import sys
import pdb
import random
from numpy import unique
from collections import defaultdict
from typing import List, Tuple, Dict
from optparse import OptionParser


def ReadJunc(junfile: str) -> Dict:
    jundict = defaultdict(list)
    fp = open(junfile, "r")

    for line in fp.readlines():
        lline = line.rstrip().split()
        
        jundict[lline[9]].append(line.rstrip().split())

    return jundict


def ComputeMultiMapping(alist):
    #alist--->[[],[],[],[]]
    if len(alist) == 1:
        return True

    clist, plist = [], []
    for aline in alist:
        if aline[0] == "chrLT":
            clist.append(aline[3])
            plist.append(int(aline[4]))
        else:
            clist.append(aline[0])
            plist.append(int(aline[1]))

    if len(unique(clist)) != 1:
        return False
    else:
        if unique(plist)[-1] - unique(plist)[0] < 20:
            return True


def SelectOneRead(plist: List[List]) -> List:
    sumr = len(plist)
    arr = defaultdict(int)

    for lline in plist:
        if lline[0] == "chrLT":
            arr[lline[4]] += 1
        else:
            arr[lline[1]] += 1

    savepos = max(arr.items(),key=lambda ads:ads[1])[0]

    for line in plist:
        if line[0] == "chrLT" and line[4] == savepos:
            return line
        if line[3] == "chrLT" and line[1] == savepos:
            return line


def UniqOneBcPosition(reads: List[List]):
    suppfrag = len(reads)
    chrList, posList = [], []

    for read in reads:
        if read[0] == "chrLT":
            chrList.append(read[3])
            posList.append(int(read[4]))
        else:
            chrList.append(read[0])
            posList.append(int(read[1]))

    if len(unique(chrList)) != 1:
        return False

    if len(unique(posList)) == 1:
        return True
    else:
        if sorted(posList)[-1] - sorted(posList)[0] < 50:
            return True
        else:
            return False


def main():
    #pdb.set_trace()
    usage = "\n\t star-junction.fiter.py [options] <Breaks.to.gene>\n\n \
    \t###########################################################\n \
    \t#\n \
    \t# Required:\n \
    \t#     --Togene | string   Junction file of breaks to gene[star-fusion.junction_breakpts_to_genes_pass.txt]\n \
    \t###########################################################"

    parser = OptionParser(usage=usage)
    parser.add_option('--Togene', action='store', dest='Gene',
            help="[required] junction file of have mapped to gene")

    (option, args) = parser.parse_args()

    if (not option.Gene):
        parser.print_help()
        sys.exit(0)

    jundict = ReadJunc(option.Gene)
    passname = os.path.join(os.path.split(option.Gene)[0],"star-fusion.multi.map.oneqname.pass.txt")
    one_qname_multi_map_pass = open(passname, "w")
    one_qname_multi_map_fail = open(os.path.join(os.path.split(option.Gene)[0],"star-fusion.multi.map.oneqname.fail.txt"), "w")

    #remove one qname have multiple position
    for qname, plist in jundict.items():
        Saved = ComputeMultiMapping(plist)
        if Saved:
            if len(plist) == 1:
                one_qname_multi_map_pass.write("%s\n" %("\t".join(plist[0])))
            else:
                l = SelectOneRead(plist)
                one_qname_multi_map_pass.write("%s\n" %("\t".join(l)))
        else:
            for line in plist:
                one_qname_multi_map_fail.write("%s\n" %("\t".join(line)))

    one_qname_multi_map_pass.close()
    one_qname_multi_map_fail.close()

    one_qname_multi_map_pass_ = open(passname, "r")
    one_bc_multi_map_pass = open(os.path.join(os.path.split(option.Gene)[0],"star-fusion.multi.map.onebc.pass.txt"),"w")
    one_bc_multi_map_fail = open(os.path.join(os.path.split(option.Gene)[0],"star-fusion.multi.map.onebc.fail.txt"),"w")

    #remove on barcode have multiple position
    bcDict = defaultdict(list)
    for line in one_qname_multi_map_pass_.readlines():
        lline = line.rstrip().split()
        qname, bc = lline[9].split("|")[0], lline[9].split("|")[1]
        bcDict[bc].append(lline)

    for bc, reads in bcDict.items():
        if UniqOneBcPosition(reads):
            one_bc_multi_map_pass.write("%s\t%d\n" %("\t".join(reads[random.randint(0,len(reads)-1)]),len(reads)))
        else:
            for read in reads:
                one_bc_multi_map_fail.write("%s\n" %("\t".join(read)))



if __name__ == '__main__':
    main()
