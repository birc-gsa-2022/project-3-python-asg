#####################
# Libraries:

import sys
import random
import math

#####################
# Functions:

from cigar import edits_to_cigar
from align import get_edits
    
# Simulate string:
def simulate_string(m):
    """Simulate a DNA sequence of length m.
    
    Example:
    simulate_string(20)
    
    """
    DNA = 'ACGT'
    nucleotides = [random.choice(DNA) for _ in range(m)]
    return ''.join(nucleotides)
    

# Naive approach
def get_suffix_array(string):
    string+='$'
    SA = sorted(range(len(string)), key=lambda i: string[i:])
    return SA

    
# Prefix-doubling approach (SIMPLE EDITION - no radix)
def SuffixArray(string):
    '''Algorithm for constructing suffix-arrays
    
    Example:
    >>> SuffixArray('trollolol')
    [9, 8, 3, 6, 4, 7, 2, 5, 1, 0]
    
    '''
    if string == '' or string == None:
        return [0]
    
    string += '$'
    rank_list = list(range(len(string)))
    index = {v: i for i, v in enumerate(sorted(set(string)))}
    string = [index[v] for v in string]
    M,j = 0,1
    while M < len(string)-1:
        tuple_list = []
        for i in range(len(string)):
            ij = i+j
            if ij < len(string):
                key = (string[i], string[ij])
            else: 
                key = (string[i], string[-1])
            tuple_list.append(key)
        j=j*2
        keys = sorted(set(tuple_list))
        ranks = dict(zip(keys, rank_list))
        string = [ranks[tuple_list[i]] for i in range(len(string))]
        M = max(string)
    SA = [None]*len(string)
    for i in range(len(string)):
        SA[string[i]] = i
    return SA


def binary_search(SA, string, pattern):
    '''Algorithm for finding exact patterns from suffix array.
    
    Example:
    string = 'trollolol'
    pattern = 'ol'
    SA = SuffixArray(string)
    binary_search(SA, string, pattern)
    [7, 2, 5]
    
    '''
    if pattern == '' or pattern == None:
        return None
    
    if string == '' or string == None:
        return None
    
    if SA == '' or SA == None:
        return None
    
    string+='$'
    SA_positions = []
    # Binary search.
    upper, lower = len(SA), 0, 
    last, match =  None, None
    j = 0
    S = 0
    max_runs = math.ceil(math.log2(len(SA)))
    while j < max_runs and S == 0:
        mid = lower + (upper - lower) // 2
        if mid == last: 
            S = 1
        count = 0
        for i in range(count, len(pattern)):
            if pattern[i] == string[SA[mid]+i]:
                count+=1
            if count == len(pattern):
                match = mid    
                SA_positions.append(SA[match])
                S = 1
                break
            elif pattern[i] > string[SA[mid]+i]:
                    lower = mid
                    break
            elif pattern[i] < string[SA[mid]+i]:
                    upper = mid
                    break
        last=mid
        j+=1 
    if match != None:
        k=1
        while match-k >= 0 and string[SA[match-k]:SA[match-k]+len(pattern)] == pattern:
            SA_positions.append(SA[match-k])
            k+=1
        k=1
        while match+k < len(string) and string[SA[match+k]:SA[match+k]+len(pattern)] == pattern:
            SA_positions.append(SA[match+k])
            k+=1  
    return SA_positions


def read_fasta():
    # load input:
    inFile = sys.argv[1]
    with open(inFile,'r') as f:
        lines = f.readlines()
    record_list = []
    header = ''
    sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if header != "":
                record_list.append([header.strip(), ''.join(sequence).strip()])
                sequence = []
            header = line[1:]
        else:
            sequence.append(line)
    record_list.append([header.strip(), ''.join(sequence).strip()])
    return record_list


def read_fastq():
    inFile = sys.argv[2]  
    with open(inFile,'r') as f:
        lines = f.readlines()
    record_list = []
    header = ''
    sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith('@'):
            if header != "":
                record_list.append([header.strip(), ''.join(sequence).strip()])
                sequence = []
            header = line[1:]
        else:
            sequence.append(line)
    record_list.append([header.strip(), ''.join(sequence).strip()])
    return record_list

#####################
# Code:
    
# if __name__ == '__main__':

#     fasta_recs = read_fasta()
#     fastq_recs = read_fastq()
    
#     for fa_rec in fasta_recs:
#         ref = fa_rec[1]
#         SA = SuffixArray(ref)
#         for fq_rec in fastq_recs:
#             read = fq_rec[1]
#             matches = binary_search(SA, ref, read)
#             for match in matches:
#                 read_name = fq_rec[0]
#                 read_seq = fq_rec[1]
#                 edits = get_edits(read_seq, fa_rec[1][match:match+len(fq_rec[1])])
#                 cigar = edits_to_cigar(edits[2])
#                 output = [read_name,fa_rec[0],str(match+1),cigar,read_seq]
#                 print('\t'.join(output))

###############################################################


string = ''
print(SuffixArray(string))
print(get_suffix_array(string))

