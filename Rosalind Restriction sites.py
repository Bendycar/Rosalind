#!/usr/bin/env python

# Solves the following problem: https://rosalind.info/problems/revp/
input = "C:/Users/Ben/Downloads/rosalind_revp.txt"
complement = {"A" : "T", "T" : "A", "G" : "C", "C" : "G", "a" : "t", "t" : "a", "g" : "c", "c" : "g", "N" : "N", "n" : "N"}

def reverse_complement(DNA: str) -> str:
    RC = ""
    for nucleotide in reversed(DNA):
        RC += complement[nucleotide]
    return RC

def check_reverse_palindrome(DNA: str) -> bool:
    if len(DNA) >= 4 and len(DNA) <= 12: #This line is unnecessary when used in main(), but is left in case this is ever used elsewhere. Problem states length must be between 4 and 12bp.
        return DNA == reverse_complement(DNA)
    return False
    
def extract_fasta(read) -> str:
    '''Takes a FASTA file and returns a string containing just the nucleotides'''
    write = ""
    with open(read, "r") as fh:
        for line in fh: 
            if line.startswith(">"): #We want to exclude the header lines
                continue
            else:
                write += (line.strip('\n')) 
    return write

FASTA = extract_fasta(input)

def main(FASTA):
    with open("C:/Users/Ben/Downloads/rosalind_revp_ans.txt", "w") as answer:
        for start in range(0, (len(FASTA))): 
            for length in range(4,13):
                end = start + length
                if end > len(FASTA):
                    break
                DNA = FASTA[start:end]
                if check_reverse_palindrome(DNA):
                    answer.write(f"{start + 1}\t{len(DNA)}\n") #Problem requested 1-based indexing

main(FASTA)






