import sys
inputfile = sys.argv[1]
fin = open(inputfile,'r')
inputsequence = ''
for line in fin:
     inputsequence = inputsequence + line.strip('\n')
     inputsequence = inputsequence.replace(' ', '')

from copy import deepcopy as dc
import random
fwd_table={'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'TAA': '.', 'TAG': '.', 'TGA': '.'}


backtable = {}
for key in fwd_table:
    val = fwd_table[key]
    if not val in backtable:
        backtable[val] = [key]
    else:
        backtable[val].append(key)
rk_switch_backtable = {'F': ['TTT', 'TTC'],
 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
 'Y': ['TAT', 'TAC'],
 'C': ['TGT', 'TGC'],
 'W': ['TGG'],
 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
 'H': ['CAT', 'CAC'],
 'Q': ['CAA', 'CAG'],
 'R': ['AAA', 'AAG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
 'I': ['ATT', 'ATC', 'ATA'],
 'M': ['ATG'],
 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
 'N': ['AAT', 'AAC'],
 'K': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG', 'AAA', 'AAG'],
 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
 'D': ['GAT', 'GAC'],
 'E': ['GAA', 'GAG'],
 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
 '.': ['TAA', 'TAG', 'TGA']}

bad_codons = ['CTA', 'ATA', 'TAG', 'CGA', 'AGA', 'AGG', 'GGA']
surface = [13, 15, 20, 23, 28, 38, 41, 44, 45, 46, 47, 50, 56, 60, 63, 70, 71, 79, 82, 86, 87, 91]
r_k_codons = ['AAA', 'AAG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
codons_not_to_mutate = []


class Triplet:
    def __init__(self,trip, verbose = False):
        self.verbose = verbose
        self.trip = trip.upper()
        assert(len(self.trip)==3) , 'Triplet not of length 3.'
        self.AA = fwd_table[self.trip]
        self.altcodons = dc(backtable[self.AA])
        self.expaltcodons = dc(rk_switch_backtable[self.AA])
        for item in self.altcodons: 
            if item in bad_codons: self.altcodons.remove(item)
        for item in self.expaltcodons: 
            if item in bad_codons: self.expaltcodons.remove(item)
        if self.trip in self.altcodons and len(self.altcodons)>1: self.altcodons.remove(self.trip)
        if self.trip in self.expaltcodons and len(self.expaltcodons)>1: self.expaltcodons.remove(self.trip)
        if self.verbose: print ('Triplet of %s initiated with amino acid being %s and alternative codons %s' %(self.trip,self.AA,self.altcodons))
    
    def to_mutate(self):
        for self.trip in codons_not_to_mutate:
            return False
        else:
            return True
        
    def random_alt(self):
        if self.to_mutate():
            random.seed()
            rannum = random.randint(0,len(self.altcodons)-1)
            return self.altcodons[rannum]
        else:
            return self.trip
    
    def random_expalt(self):
        if self.to_mutate():
            random.seed()
            rannum = random.randint(0,len(self.expaltcodons)-1)
            return self.expaltcodons[rannum]
        else:
            return self.trip
       
    
def split_seq_to_Triplets(sequence, verbose = False):
    assert(len(sequence)%3 ==0), 'Sequence length not divisible by 3. Not gonna work with codon triplets.'
    list_of_triplet_objects = []
    for i in range(0,len(sequence),3):
        list_of_triplet_objects.append(Triplet(sequence[i:i+3],verbose))
    return list_of_triplet_objects

def new_random_sequence(sequence, verbose = False):
    newseq = ""
    triplets = split_seq_to_Triplets(sequence, verbose)
    for i in range(0, len(triplets)):
        if i+2 in surface:
            newseq = newseq + triplets[i].random_expalt()
        else:
            newseq = newseq + triplets[i].random_alt()
    if verbose: print ('Writing old sequence %s to new sequence %s' %(sequence,newseq))
    return newseq

def simple_compare(seq1,seq2):
    assert(len(seq1)==len(seq2)), 'Sequence 1 not of same length as sequence 2. Obviously makes no sense like this.'
    hits = 0
    for i,c in enumerate(seq2):
        if c.lower() == seq1[i].lower():
            hits += 1
    percent = round(hits/len(seq1)*100,2)
    return hits, percent

def ortho_sequences(sequence):
    scorelist = []
    newseqlist = []
    for num in range(0,20):
        newseq = new_random_sequence(sequence)
        score = simple_compare(newseq, sequence)[1]
        scorelist.append(score)
        newseqlist.append(newseq)
    return scorelist, newseqlist, newseqlist[scorelist.index(min(scorelist))]
    
newseq = ortho_sequences(inputsequence)  

"""for num in range(0,10):
    newseq = ortho_sequences(inputsequence)
    with open(inputfile.split('.')[-2]+'_new_'+str(num)+'.fasta','w') as fout:
       fout.write(newseq[0])
       fout.write(newseq[1])
       fout.write(newseq[2])"""

"""Useful stuff"""
'''raw = 'gctgaccgtatcgaactgcgcggcctgaccgtgca tggtcgccac ggggtctacg accacgagcg tgtggccggg cagcgctttg tcatcgatgt caccgtgtgg attgacctgg ccgaggccgc caacagcgac gacttggccg acacctatga ctacgtgcgc ctggcttcgc gcgcggccga gatcgtcgcc ggacccccgc gcaagctgat cgaaacggtc ggggccgaga tcgctgatca cgtgatggac gaccagcgtg tgcatgccgt tgaggtggcg gtacacaagc cgcaggcgcc cattccgcag acgttcgacg atgtggcggt ggtgatccga cgctcacggc gctaa'
gene = ''
for c in raw:
    if not c == ' ':
        gene = gene + c
new = new_random_sequence(gene[:-3])
'''
