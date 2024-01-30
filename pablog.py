import sys
from traceback import print_exc
import pysam
import HTSeq as ht

from utils import utilsAlignment 
from utils import otherUtils
from utils import utilsPrimer3


def main(argv):

   if(len(argv) < 4) or (len(argv) > 5):
        print("\nRun pablog.py script using python:\n")
        print("   python pablog.py FILE.bam Ref_Gene.fa output_result_file.txt [optional]\n")
        print(" where:\n - FILE.bam is the result file from alignment,\n")
        print(" - Ref_Gene.fa is the reference genome used in the alignment,\n")
        print(" - output_result_file.txt filename where results are written in,\n")
        print(" - [optional] args are: \n      - size [INT] (=60 by default) which is a the minimun length we consider to filter out cigar operations in read analysis.\n")
        return 
   elif(len(argv) == 4):
         
         
         try:
             #BAM file parsing
             aln_file = ht.BAM_Reader(argv[1])
         
         except Exception as err:
             
             print(err)
             return
         
         print("Analyzing alignment reads... \n")
         regions = utilsAlignment.analyzeAlignment(aln_file)   #Alignment Analysis step
         print("DONE!\n\n")

         if len(regions) > 0 :
            
            print("Analyzing regions' coverage ... \n")
            sorted_regions = utilsAlignment.analyzeCoverage(argv[1],regions)  #Coverage Analysys step 
            print("DONE!\n\n")

            try:

               print("Generating Consensus sequence from bam file ...\n")
               consensus_sequence = utilsPrimer3.generateConsensus(argv[1],argv[2]) #Consensus Sequence step
               print("DONE!\n\n")
               
               
               print("Primers design ...\n")
               utilsPrimer3.p3Design(consensus_sequence, regions, argv[3])    #Primers Design step 
               print("DONE!\n")

               print("PABLOG analysis completed, results available in %s" %argv[3])

               return
            
            except Exception as err:
               print(err)
               return
         
         else: 

            print("No good regions were found from the alignment.\n"
                  "No primer was designed.\n"
                  "Check if alignment is consistent or try using less restricting parameters.\n")
            
            return

   else:
         
         #Same pipeline as above, but with optional parameter region min-length provided 

         try:
             
             aln_file = ht.BAM_Reader(argv[1])
         
         except Exception as err:
             
             print(err)
             return
         
         print("Analyzing alignment reads... \n")      
         regions = utilsAlignment.analyzeAlignment(aln_file, _size=argv[4])
         print("DONE!\n\n")

         if len(regions) > 0:

            print("Analyzing regions' coverage ... \n")
            sorted_regions = utilsAlignment.analyzeCoverage(argv[1],regions)
            print("DONE!\n\n")

            try:
               
               print("Generating Consensus sequence from bam file ...\n")
               consensus_sequence = utilsPrimer3.generateConsensus(argv[1],argv[2]) #Consensus Sequence step
               print("DONE!\n\n")

               print("Primers design ...\n")
               utilsPrimer3.p3Design(consensus_sequence, regions, argv[3])    #Primers Design step 
               print("DONE!\n\n")

               print("PABLOG analysis completed, results available in %s" %argv[3])

               return
            
            except Exception as err:

               print(err)
               return
         else: 

            print("No good regions were found from the alignment.\n"
                  "No primer was designed.\n"
                  "Check if alignment is consistent or try using less restricting parameters.\n")
            return


if __name__ == "__main__":
#Main function
   try:
      arg = sys.argv
      main(sys.argv)

   except Exception as err:
      print("Error: ", err)
      sys.exit(-1)
