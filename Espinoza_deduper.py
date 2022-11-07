#!/usr/bin/env python3.10

import argparse
import re
import gzip


#Make a Note that the user will require a sorted and paired SAM file prior to using this script 
def get_args():
    parser = argparse.ArgumentParser(description="A program for PCR deduplication to obtain accurate RNA seqs")
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file", required=True,)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to sorted sam file", required=True,)
    parser.add_argument("-u", "--umi", help="Set of known UMIs", required=True, type=str)
    return parser.parse_args()

args=get_args()

#---------- Dictionaries to count each condition ------------
#first store all recods inside a dictionary and then, any repeating records add to a deduplicate_count
#unique_reads:{(chromosome, adj_pos, flag, umi)}
unique_reads = {}

##---------- UMI Set ------------
# save known UMIs inside of a set to check if UMIS in SAM file match them and if not = trash
umi_set = set()

##---------- Counters ------------
wrong_umis = 0
header_count = 0
unique_count = 0
deduplicates_count = 0

#This is to open STL96.txt file anf store each UMI in a set
with open(args.umi, 'r') as known_umis:
    for line in known_umis:
        umi_set.add(line.strip())
    

#---------- Functions ------------
#Strand Directionality 
def strand_direc(flag:int)->bool:
    Reverse = True
    if ((flag & 16) == 16):
        return Reverse
    else:
        Reverse = False
        return Reverse

#---------- Finding a cigar string pattern and storing each into a set ------------
def cigar_search(cigar:str) -> list[str]:
    cigar_pattern="[0-9]+[MIDNSHPX]+"
    cigar_list=re.findall(cigar_pattern, cigar)
    return cigar_list


#---------- Correcting for start position ------------
# cigar [str] | start_position of strand [str] | strand direction [bool] 
def adjusted_position(cigar_list:list[str], position:int, Reverse:bool)-> int:  # type: ignore
   
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
        elif letter == "I":
            insertions+= number
        elif letter == "D":
            deletions += number
        elif letter == "N":
            splice += number
        elif letter == "S": 
            if letter in cigar_list[0]:
                left_soft_clipped += number
            elif letter in cigar_list[-1]:
                right_soft_clipped += number

       
#----- forward and reverse adjusted position:
        #Forward = False
 
    if flag == False:
        adjusted_position=(int(position) - left_soft_clipped)
        #Reverse = True
    elif flag == True:
        adjusted_position=(int(position) + matches + deletions + splice + right_soft_clipped)
    return adjusted_position


#---------- For Loop ------------
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
           
            Reverse=strand_direc(flag)
            adj_pos=adjusted_position(cigar_list, position, Reverse) 
            my_key=(umi, flag, chromosome, adj_pos)
        
            if my_key not in unique_reads:
                unique_count += 1
                unique_reads[my_key]=1
                outfile.write(line)
            else:
                deduplicates_count+= 1
            

print("Number of Header Lines: ", header_count)
print("Number of Unique Reads: ", unique_count)
print("Number of Deduplicates: ", deduplicates_count)
print("Number of Wrong Umis: ", wrong_umis)




