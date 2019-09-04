###################################################################
# File Name: pipe.sh
# Author: yys
# mail: shayy0919@163.com
# Created Time: 2019年09月03日 星期二 10时36分53秒
###################################################################
#!/bin/bash

python3 star-junction.Match_length.py \
    --Amplicon /home/yys/Project/Lilab/Lentivirus/Data/SecondPrimData/LG.primer \
    --Chimeric /home/yys/Project/Lilab/Software/Result/Chimeric.out.junction \
    --Alignbam /home/yys/Project/Lilab/Software/Result/STAR/Aligned.out.bam 

python3 star-junction.map_chimeric_reads_to_genes.py \
    --JuncM /home/yys/Project/Lilab/Software/Result/std.Chimeric.out.junction.Mathch_length \
    --Gene /home/yys/Project/Lilab/Lentivirus/DataBase/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans \
    --Exon /home/yys/Project/Lilab/Lentivirus/DataBase/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/ref_annot.gtf.mini.sortu.tmp

python3  star-junction.fiter.py \
    --Togene /home/yys/Project/Lilab/Software/Result/star-fusion.junction_breakpts_to_genes_pass.txt
