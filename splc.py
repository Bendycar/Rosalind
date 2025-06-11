#!/usr/bin/env python
# https://rosalind.info/problems/splc/

import argparse

codons = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',    # Serine
    'TTC': 'F', 'TTT': 'F',    # Phenylalanine
    'TTA': 'L', 'TTG': 'L',    # Leucine
    'TAC': 'Y', 'TAT': 'Y',    # Tyrosine
    'TAA': '*', 'TAG': '*',    # Stop
    'TGC': 'C', 'TGT': 'C',    # Cysteine
    'TGA': '*',    # Stop
    'TGG': 'W',    # Tryptophan
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',    # Leucine
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',    # Proline
    'CAC': 'H', 'CAT': 'H',    # Histidine
    'CAA': 'Q', 'CAG': 'Q',    # Glutamine
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',    # Arginine
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',    # Isoleucine
    'ATG': 'M',    # Methionine
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',    # Threonine
    'AAC': 'N', 'AAT': 'N',    # Asparagine
    'AAA': 'K', 'AAG': 'K',    # Lysine
    'AGC': 'S', 'AGT': 'S',    # Serine
    'AGA': 'R', 'AGG': 'R',    # Arginine
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',    # Valine
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',    # Alanine
    'GAC': 'D', 'GAT': 'D',    # Aspartic Acid
    'GAA': 'E', 'GAG': 'E',    # Glutamic Acid
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'     # Glycine
}

def get_args():
    parser = argparse.ArgumentParser(description="Rosalind problem 'splc' -- deletes introns and translates exons from provided FASTA file. File is formatted with the RNA string in the first entry and introns in all other entries.")
    parser.add_argument("-i", "--input", help="FASTA file from Rosalind website", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output file to store protein sequence rather than printing to stdout, if desired", type=str, required=False)
    return parser.parse_args()

args = get_args()
input_file = args.input
output_file = args.output

def parse_fasta(fasta_file):
    '''Takes in FASTA file as input, returns tuple containing the RNA sequence and list of introns'''

    with open(input_file, 'r') as fh:
        introns = []
        fh.readline() # Getting rid of first sequence header 
        sequence = '' 
        first = True # Necessary to parse potentially multi-line sequence in first entry
        while first:
            line = fh.readline().strip()
            if line.startswith('>'): # Meaning we have reached the first intron entry
                first = False
                break
            else:
                sequence += line

        for line in fh:
            if line.startswith('>'):# Sequence header
                continue
            else:
                introns.append(line.strip())
    
    for intron in introns:
        sequence = sequence.replace(intron,'') # Delete intron by replacing with empty string

    return(sequence, introns)

def translate(sequence, introns):
    '''Given full RNA sequence and list of introns, removes introns from sequence and translates remaining exons'''

    protein_sequence = ''

    for i in range(0, len(sequence) - 3, 3): 
        codon = sequence[i:i+3]
        protein_sequence += codons[codon]

    return(protein_sequence)

if __name__ == "__main__":
    sequence, introns = parse_fasta(input_file)
    protein_sequence = translate(sequence, introns)

    if output_file:
        with open(output_file, "w") as ofh:
            ofh.write(protein_sequence)
    else:
        print(protein_sequence)


