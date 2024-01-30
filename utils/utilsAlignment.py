import pysam
import pandas as pd
import numpy as np
import tqdm as tq
import time
from traceback import print_exc

from .classRegion import Region
from .classRegion import addUniqueRegion



accepted_types = ['N', 'H', 'S']

#analyzeRead takes an aligned read from alignment file 
def analyzeRead(read_alnd, _size):
  
    #filter out cigar ops sized lower than given _size (default 60)

    filtered_ops = []
    for op in read_alnd.cigar:
        if(op.size > _size):
          filtered_ops.append(op)

    number_of_ops = len(filtered_ops)

    #the result is a list of intervals, deduced from cigar ops 
    good_cigar_ops = [] 

    for i,op in enumerate(filtered_ops):
        if(i != 0 and i < number_of_ops - 1):
          prev_op = filtered_ops[i-1]
          next_op = filtered_ops[i+1]

          if(op.type in accepted_types) and (next_op.type == prev_op.type == 'M'):
              good_cigar_ops.append(op)
      
    if len(good_cigar_ops) > 0:
        return [read_alnd.read.name, good_cigar_ops, read_alnd.aQual]
    
    return []


def analyzeAlignment(aln_file, _size=60):
#Implements Alignment Analysis step
    
    #list of unique regions
    unique_regions = []

    try:
      
      for line in tq.tqdm(aln_file):
        
        if(line.aligned):
            
            result = analyzeRead(line, int(_size))

            if(result):
                for element in result[1]:
                    start_ = element.ref_iv.start
                    end_ = element.ref_iv.end
                    name_ = element.ref_iv.chrom
                    region = Region(s=start_, e=end_, name=name_) 
                    addUniqueRegion(unique_regions, region)
        else:
           break
    except Exception as err:
       print_exc()
    
    return unique_regions

def analyzeCoverage(bamfilename, regions):
  #Implements the coverage analysis step
  #Goodness Score is computed for each candidate region found, using coverage variation formula

  try:
      gene_coverage = pysam.depth(bamfilename, "-a")

  except Exception as err:
      print(err)
      return
  
  data_frame = pd.DataFrame([x.split('\t') for x in gene_coverage.split('\n')])
  data_frame.columns = ['chr', 'index', 'coverage']

  data_frame[["index","coverage"]] = data_frame[["index","coverage"]].apply(pd.to_numeric)

  for region in tq.tqdm(regions):
      
      local_start_index = region.start - 1 
      local_end_index = region.end + 1 
      interval_dim = local_end_index - local_start_index

      df_slice = data_frame.iloc[local_start_index:local_end_index]
      c_m1 = df_slice.iloc[0,2]
      c_m2 = df_slice.iloc[interval_dim - 1,2]
      c_max = df_slice.iloc[1:interval_dim - 2].max()['coverage']

      mean_coverage = (c_m1 + c_m2) / 2
      variation = mean_coverage - c_max
      variation_percentage = (variation/mean_coverage) * 100
      
      if not (variation_percentage > 0.0):
      	region.set_confidence(0.0)
      else:	
      	region.set_confidence(variation_percentage)
  
  sorted_regions = sorted(regions, key = lambda x: x.confidence , reverse = True)
    
  return sorted_regions
   


