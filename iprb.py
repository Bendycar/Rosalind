#!/usr/bin/env python
# https://rosalind.info/problems/iprb/

import random
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Rosalind Mendelian inheritance simulator. Takes in counts representing homozygous dominant, heterozygous, and homozygous recessive individuals in a population"
)
    parser.add_argument("-k", help="Homozygous dominant", type=int, required=True)
    parser.add_argument("-m", help="Heterozygous", type=int, required=True)
    parser.add_argument("-n", help="Homozygous recessive", type=int, required=True)
    return parser.parse_args()

args = get_args()
k = args.k 
m = args.m 
n = args.n 

def simulation(k,m,n):
    
    results_A = 0 # Dominant phenotype
    results_a = 0 # Recessive phenotype

    for i in range(1000000): # We want to simulate this probability with a very high sample size, rather than explicitly calculate it. TODO: Add explicit calculation
        
        mates = random.sample(['k', 'm', 'n'], counts = [k, m, n], k=2) # Select 2 individuls from the population according to their relative proportions
        
        if mates == ['n', 'n']: # Two recessive parents
            results_a += 1
            
        elif mates == ['m', 'n'] or mates == ['n', 'm']:# Heterozygous and recessive parents
            
            offspring_mn = random.sample(['A', 'a'], k = 1) # Even probability of offspring phenotype, according to Punnett square
            if offspring_mn == ['a']:
                    results_a += 1
            else: 
                    results_A += 1
                
        elif mates == ['m', 'm']: # Two heterozygous parents
                
            offspring_mm = random.sample(['A', 'a'], counts = [75, 25], k = 1) # Punnett Square probability again
            if offspring_mm == ['a']:
                    results_a += 1
            else: 
                    results_A += 1
    
        else: # Two homozygous dominant parents will always return dominant phenotype
                    results_A += 1

    return results_A / (results_A + results_a)

print(simulation(k,m,n))
