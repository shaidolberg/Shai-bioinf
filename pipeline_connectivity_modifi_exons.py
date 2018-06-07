import os
import re

inFile_mod = open(r"/home/shaidulberg/chipseq/Modifications/mod_list.txt","r")
line = inFile_mod.readline()
d={}
while line:
    if re.search('pval',line.rstrip()):
        ls = line.rstrip().split()
        line = inFile_mod.readline()
        #print('line\n')
        #os.system(line)
        m = re.search('/@@download/(.*)\.bigWig$',line.rstrip())
        d[ls[1]] = m.group(1)
    line = inFile_mod.readline()

#os.system('module load bwtool\n')

for modification in d:
 #   print('bwtool extract bed -tabs /home/shaidulberg/chipseq/input/for_bwtool_exons/exons_up_for_bwtool.txt  /home/shaidulberg/chipseq/bigWigs/'+d[modification]+'.bigWig   /home/shaidulberg/chipseq/BigWig_extract/exons/'+modification+'.pval_exon.table_up.txt')
  #  os.system('bwtool extract bed -tabs /home/shaidulberg/chipseq/input/for_bwtool_exons/exons_up_for_bwtool.txt  /home/shaidulberg/chipseq/bigWigs/'+d[modification]+'.bigWig   /home/shaidulberg/chipseq/BigWig_extract/exons/'+modification+'.pval_exon.table_up.txt')
    
   # print('bwtool extract bed -tabs /home/shaidulberg/chipseq/input/for_bwtool_exons/exons_down_for_bwtool.txt  /home/shaidulberg/chipseq/bigWigs/'+d[modification]+'.bigWig  /home/shaidulberg/chipseq/BigWig_extract/exons/'+modification+'.pval_exon.table_down.txt')
    #os.system('bwtool extract bed -tabs /home/shaidulberg/chipseq/input/for_bwtool_exons/exons_down_for_bwtool.txt  /home/shaidulberg/chipseq/bigWigs/'+d[modification]+'.bigWig   /home/shaidulberg/chipseq/BigWig_extract/exons/'+modification+'.pval_exon.table_down.txt')

#os.system('module load  R/R-3.4.3\nmodule load gcc/gcc630\n')
    print('Rscript /home/shaidulberg/chipseq/Modifications/line_plot_by_connectivity_unix.R /home/shaidulberg/chipseq/BigWig_extract/exons_clean/'+modification+'.pval_exon.table_up_clean.txt   /home/shaidulberg/chipseq/BigWig_extract/exons_clean/'+modification+'.pval_exon.table_down_clean.txt')   
    os.system('Rscript /home/shaidulberg/chipseq/Modifications/line_plot_by_connectivity_unix.R /home/shaidulberg/chipseq/BigWig_extract/exons_clean/'+modification+'.pval_exon.table_up_clean.txt   /home/shaidulberg/chipseq/BigWig_extract/exons_clean/'+modification+'.pval_exon.table_down_clean.txt')    
    
                    
inFile_mod.close()       
