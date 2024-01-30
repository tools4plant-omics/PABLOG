from .classRegion import Region 


def count_occurrencies(string):
#Utility for consensus sequence 
    counters = [0,0,0,0,0]
    ret = ['ref','A','C','G','T']
    #          .,A,C,G,T                
    for char in string:
        if char =='.' or char==',' or char=='<' or char=='>':
            counters[0] +=1
        elif char =='A' or char=='a':
            counters[1]+=1
        elif char =='C' or char=='c':
            counters[2]+=1
        elif char =='G' or char == 'g':
            counters[3]+=1
        elif char == 'T' or char == 't':
            counters[4]+=1
    index = counters.index(max(counters))
    return ret[index]

def get_base_from_line(line,stringa):
#Utility for consensus sequence
    splitter = line.split("\t")
    base = count_occurrencies(splitter[4])

    if base == 'ref':
        stringa+=(splitter[2])
    else:
        stringa+=base
        
    return stringa

def writeRegions(list_of_regions, outpute_filename):
    f = open(outpute_filename, "w")
    for region in list_of_regions:
        f.write(str(region) + "\n")
    f.close()

def check_overlapping_regions(regions):
    regions_ret = []
    i = 0
    end = regions[i][1]
    while(i >= 0 and i<(len(regions)-1)):
        if(regions[i][0] < regions[i+1][1]):
            i+=1
        else:
            start = regions[i][0]
            regions_ret.append((start,end))
            i+=1
            end = regions[i][1]
    regions_ret.append((regions[i][0],end))

    return regions_ret

def check_overlapping_regions_v2(regions):
#In the process of removing introns, overlapping regions has to be managed correctly
#this function solves the problem 

    i = 0
    regions_ret = []
    shifted_indices = []

    shift_counter = 0

    while(i<(len(regions)-1)):
        #if(regions[i][0] == regions[i+1][1]): #regions with same start index -> we take the longer one 
            #if(regions[i][1] < regions[i+1][0]):
            #i+=1
        #else:
        shifted_indices.append(regions[i][0] - shift_counter)
        regions_ret.append((regions[i][0] - shift_counter,regions[i][1] - shift_counter))
        shift_counter += (regions[i][1] - regions[i][0])
        i+=1

    shifted_indices.append(regions[i][0] - shift_counter)
    regions_ret.append((regions[i][0] - shift_counter,regions[i][1] - shift_counter))
    return regions_ret,shifted_indices


def clean_consensus(sequence):
#CLEAN CONSENSUS STRING
#in the first version consensus sequence was retrieved by a file and needed cleaning from \n \r
#now we generate it directly into a clean string, no needed anymore

    for i in range(len(sequence)):
        if sequence[i] == '\n' or sequence[i] == '\r':
            break

    tmp_seq=sequence[i:]
    seq_ret = tmp_seq.replace('\n', '')

    return seq_ret


def sequence2primer3(sequence,bed_regions):
#EDIT CONSENSUS SEQUENCE FOR PRIMER3 INPUT
#we want to cut away substrings corresponding to regions 
#(they are introns, we don't want them for RT-qPCR experiments... )

    seq = sequence
    regions=[]

    for item in bed_regions:
        regions.append((item.start,item.end))
        
    regions_sorted = sorted(regions, key=lambda x: (x[0],x[1]))
    regions_sorted,indices = check_overlapping_regions_v2(regions_sorted)
    
    for item in regions_sorted:
        seq = seq.replace(seq[item[0]:item[1]],'')

    return seq,indices

def p3Args_fromfile(filename):
#Primer3 default settings are parsed from a static file inside PABLOG folder (primer3_settings.txt)
#User can edit that file for specific values according to needings

    p3args = {}
    with open(filename) as f:
        for line in f:
            (key, val) = line.split()
            p3args[str(key)] = val

    for key,value in p3args.items():
        if not isinstance(value,int):
            if '.' in value:
                try:
                    p3args[key] = float(value)
                except ValueError:
                    pass
            else:
                try:
                    p3args[key] = int(value)
                except ValueError:
                    pass

    #single paramenter setting DIRTY WAY!
    p3args['PRIMER_PRODUCT_SIZE_RANGE'] = [80, 150]

    return p3args

def print_info_p3(dict_data, fileIO):
#Takes Primer3 results (dict format) and print into results file 
    
    data_keys = ['TM', 'GC_PERCENT','SELF_ANY_TH', 'SELF_END_TH','HAIRPIN_TH', 'SEQUENCE']
    data_left = dict_data['PRIMER_LEFT'][0]
    data_right = dict_data['PRIMER_RIGHT'][0]

    to_print = ''
    to_print+='OLIGO \t\t\t'+ 'start'+'\t\t'+'len'+'\t\t'+'tm'+'\t\t'+'gc%'+'\t\t'+'any_th'+'\t\t'+'3\'_th'+'\t\t'+'hairpin'+'\t\t'+'seq'+"\n"
    to_print+=( "LEFT PRIMER \t\t"+str(data_left['COORDS'][0]) + '\t\t'+str(data_left['COORDS'][1])+'\t\t')

    for key in data_keys:
        if (isinstance(data_left[key],float)):
            to_print +=(f'{data_left[key]:.2f}'+'\t\t')
        else:
            to_print+=(str(data_left[key])+'\t\t')
    to_print+='\n'

    to_print+=( "RIGHT PRIMER \t\t"+str(data_right['COORDS'][0]) + '\t\t'+str(data_right['COORDS'][1])+'\t\t')

    for key in data_keys:
        if (isinstance(data_right[key],float)):
            to_print+=(f'{data_right[key]:.2f}'+'\t\t')
        else:
            to_print+=(str(data_right[key])+'\t\t')
    to_print+='\n'

    return to_print
    
    
