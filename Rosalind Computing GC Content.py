# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 18:43:38 2024

@author: bendc
"""

def GC(string):
    GC = 0
    for i in string:
        if i == 'C' or i == 'G':
            GC += 1
    
    return (GC / len(string)) * 100
    

def extractDNA(sequences):
    with open(sequences, 'r') as sequence:
        dna = {}
        
        for line in sequence:
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
                dna[key] = ''
            else:
                dna[key] += line
        
        return dna

def maxGC(sequences):
    maxGCcount = 0
    dna = extractDNA(sequences)
    maxKey = ''
    
    for key in dna:
        if GC(dna[key]) > (maxGCcount):
            maxGCcount = GC(dna[key])
            maxKey = key
        
    return maxKey, maxGCcount
    
print(maxGC(r"C:\Users\bendc\Downloads\rosalind_gc (2).txt"))
