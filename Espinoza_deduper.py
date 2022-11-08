#!/usr/bin/env python3.10

'''
A program for PCR deduplication to obtain accurate RNA seqs.
User needs a sorted and paired SAM file prior to using this script 
'''

import argparse
import re
import gzip

def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    # parser = argparse.ArgumentParser(description='''A program for PCR deduplication to obtain accurate RNA seqs.
    # User needs a sorted and paired SAM file prior to using this script ''')
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file", required=True,)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to sorted sam file", required=True,)
    parser.add_argument("-u", "--umi", help="Set of known UMIs", required=True, type=str)
    return parser.parse_args()

args=get_args()

##---------- Dictionaries ------------
#First store all recods inside a dictionary and then, any repeating records add to a deduplicate_count
#Format: unique_reads:{(chromosome, adj_pos, flag, umi)}
unique_reads = {}

##---------- UMI Set ------------
#Save known UMIs inside of a set to check if UMIS in SAM file match them and if not = discard
umi_set = set()

##---------- Counters ------------
wrong_umis = 0
header_count = 0
unique_count = 0
deduplicates_count = 0

##---------- Opening STL96.txt ------------
#Opening STL96.txt file and storing each known UMI in umi_set()
with open(args.umi, 'r') as known_umis:
    for line in known_umis:
        umi_set.add(line.strip())
    

##---------- Functions ------------

#>>>>>>>>>> Strand Directionality <<<<<<<<<<
def strand_direc(flag:int)->bool:
    '''
    This function is used to parse the  SAM bitwise flag 16. It will check if a read is reverse complemented or not.
    It will return TRUE (set to 16) if read is in the Reverse Direction and return FALSE if read is in the Forward direction.
    '''
    bool = True
    if ((flag & 16) == 16):
        return bool
    else:
        bool = False
        return bool

# >>>>>>>>>> Searching for a CIGAR String <<<<<<<<<<
#Finding a cigar string pattern and storing each into a set
def cigar_search(cigar:str) -> list[str]:
    ''' 
    This function will search for a CIGAR string in the given SAM file and 
    add each pattern to a list called cigar_list.
    '''
    cigar_pattern="[0-9]+[MIDNSHPX]+"
    cigar_list=re.findall(cigar_pattern, cigar)
    return cigar_list


#>>>>>>>>>> Adjusting for Correct Position <<<<<<<<<<
#Format: cigar [str] | start_position of strand [str] | strand direction [bool] 
def adjusted_position(cigar_list:list[str], position:int, Reverse:bool)-> int:  
    '''
    This function uses the CIGAR pattern and adjust for the correct left most bassed start position. 
    This function takes into account cigar pattern, start position given by SAM file and strand direction.
    Insertions are ignored here.
    '''
    matches = 0
    insertions = 0
    deletions = 0
    splice = 0
    left_soft_clipped = 0
    right_soft_clipped = 0
    
    adjusted_position = 0
    for item in cigar_list:
        letter = item[-1]
       
        number = int(item[0:-1])
        if letter == "M":
            matches += number
        elif letter == "D":
            deletions += number
        elif letter == "N":
            splice += number
        elif letter == "S": 
            if letter in cigar_list[0]:
                left_soft_clipped += number
            elif letter in cigar_list[-1]:
                right_soft_clipped += number

#Forward = False
    if bool == False:
        adjusted_position=(int(position) - left_soft_clipped)

#Reverse = True
    elif bool == True:
        adjusted_position=(int(position) + matches + deletions + splice + right_soft_clipped)
    return adjusted_position


#---------- For Loop: Writing unique reads to output file ------------
with open(args.file, 'r') as fh, open(args.outfile,'w') as outfile:
    for line in fh:
        if line.startswith("@"):
            header_count += 1
            outfile.write(line) 
              
        else:
            line_split = line.split('\t') 
            readname = line_split[0]
            umi = readname.split(':')[7]
            if umi not in umi_set:
                wrong_umis += 1
                continue
    
            flag=int(line_split[1])
            chromosome=line_split[2]
            position=int(line_split[3])
            cigar=line_split[5]
            cigar_list=cigar_search(cigar)
           
            bool=strand_direc(flag)
            adj_pos=adjusted_position(cigar_list, position, bool) 
            my_key=(umi, flag, chromosome, adj_pos)
        
            if my_key not in unique_reads:
                unique_count += 1
                unique_reads[my_key]=1
                outfile.write(line)
            else:
                deduplicates_count+= 1

##---------- Print statements ------------       
print("Number of Header Lines: ", header_count)
print("Number of Unique Reads: ", unique_count)
print("Number of Deduplicates: ", deduplicates_count)
print("Number of Wrong Umis: ", wrong_umis)





