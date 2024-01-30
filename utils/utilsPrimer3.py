#This file contains functions/utilities for Primer3 step 

import os
import primer3 as p3
import pysam
from pysam import bcftools
from utils import otherUtils

from typing import (
    Any,
    Dict,
)

def basic_primer_design(seq_id, p3seq, p3args, index) -> Dict[str, Any]:
#Primer design function that wraps primer3-py call
#Return primers results as dict
    
   seq_args = {
      'SEQUENCE_OVERLAP_JUNCTION_LIST': index,
      'SEQUENCE_ID': seq_id,
      'SEQUENCE_TEMPLATE': p3seq,
   }

   design_result_dict = p3.design_primers(
      seq_args=seq_args,
      global_args=p3args,
   )

   return design_result_dict


#GENERATE CONSENSUS SEQUENCE FROM .BAM FILE
def generateConsensus(bamfile,ref):
#Implements Consensus Sequence generation step

    mpileup = pysam.mpileup('-f',ref,'-a',bamfile)
    
    with open("sam_mpileup",'w') as file:
        file.write(mpileup)
    file.close()
    
    pileup_file = open("sam_mpileup",'r')
    
    cons_seq = ''
    while True:
        line = pileup_file.readline()
        if not line:
            break
        cons_seq = otherUtils.get_base_from_line(line,cons_seq)

    
    pileup_file.close()
    os.remove("sam_mpileup")
    
    filename = bamfile[:-4] + '_consensus_sequence.fa'
    header = '>' + bamfile[:-4] + '_consensus_sequence' + '\n'
    
    with open(filename,'w') as file:
        file.write(header)
        file.write(cons_seq)
    file.close()
    
    return cons_seq


def p3Design(cons_sequence, regions, output_filename='pablog_results.txt'):
#Primer3 Design step using primer3 utility
#Sequence template input is generated and parameters parsed from settings file
#primer3-py wrap is called and results returned  
#Results are printed out in PABLOG result's text file     

   p3seq,indices = otherUtils.sequence2primer3(cons_sequence, regions)
   p3args = otherUtils.p3Args_fromfile("primer3_settings.txt")
  
   output_file = open(output_filename,'w')
   output_file.write('PABLOG ANALYSIS RESULTS...\n\n\n')
   output_file.write("Total regions found for primer designing: "+ str(len(indices)) + '\n\n')
   
   output_file.write('RefName \t'+ 'start'+'\t\t'+'end'+'\t\t'+'goodness score'+"\n")
   for region in regions:
      output_file.write(str(region) +'\n')

   output_file.write("\n\n\nFollows primers designed by Primer3 for each region found (order by goodness score): \n\n")
   
   list_of_strings_prints = []
   for i,region in enumerate(regions):
      primer = basic_primer_design(region.refname, p3seq, p3args, indices[i])
      to_print = ''
      if(primer['PRIMER_LEFT'] != []):

         to_print = ''
         to_print+=(str(region) +'\n\n')
         ret = otherUtils.print_info_p3(primer,output_file)
         to_print+=ret
         to_print+=("------------------------------------------------------------------------------\n")
         to_print+=("------------------------------------------------------------------------------\n")
         
      else:

         to_print+=(str(region) +'\n\n')
         to_print+=("NO PRIMERS DESIGNED IN THIS REGION (try different primer3 settings)\n")
         to_print+=("------------------------------------------------------------------------------\n")
         to_print+=("------------------------------------------------------------------------------\n")

      list_of_strings_prints.append([to_print,region.confidence])
   sorted_print = sorted(list_of_strings_prints, key=lambda x:x[1], reverse=True)

   for string in sorted_print:
      output_file.write(string[0])

   output_file.write('\n\n')
   output_file.write("INPUT SEQUENCE USED FOR PRIMER3:\n")
   output_file.write(p3seq)
   
   output_file.close()














