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
        return None
    string += '$'
    index = {v: i for i, v in enumerate(sorted(set(string)))}
    string = [index[v] for v in string]
    rank_list = list(range(len(string)))
    SA = [None]*len(string)
    tuple_list = SA[:]
    M,j = 0,1
    while M < len(string)-1:
        for i in range(len(string)):
            i_j = i+j
            if i_j < len(string):
                key = (string[i], string[i_j])
            else: 
                key = (string[i], string[-1])
            tuple_list[i] = key
        j*=2
        keys = sorted(set(tuple_list))
        ranks = dict(zip(keys, rank_list))
        string = [ranks[tuple_list[i]] for i in range(len(string))]
        M = max(string)
    for i in rank_list:
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
    
    # Binary search:
    hi, lo = len(SA), 0, 
    match = None
    B = 0
    while lo < hi and B == 0:
        mid = (lo + hi) // 2
        count = 0
        for i in range(count, len(pattern)):
            if pattern[i] == string[SA[mid]+i]:
                count+=1
            if count == len(pattern):
                match = mid    
                SA_positions.append(SA[match])
                B=1
            elif pattern[i] > string[SA[mid]+i]:
                lo = mid + 1
                break
            elif pattern[i] < string[SA[mid]+i]:
                hi = mid
                break
    # Scan up/down from match pos:
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
    
if __name__ == '__main__':

    fasta_recs = read_fasta()
    fastq_recs = read_fastq()
    
    for fa_rec in fasta_recs:
        ref = fa_rec[1]
        SA = SuffixArray(ref)
        for fq_rec in fastq_recs:
            read = fq_rec[1]
            matches = binary_search(SA, ref, read)
            for match in matches:
                read_name = fq_rec[0]
                read_seq = fq_rec[1]
                edits = get_edits(read_seq, fa_rec[1][match:match+len(fq_rec[1])])
                cigar = edits_to_cigar(edits[2])
                output = [read_name,fa_rec[0],str(match+1),cigar,read_seq]
                print('\t'.join(output))

###############################################################



