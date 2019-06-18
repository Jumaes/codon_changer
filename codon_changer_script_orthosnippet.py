from copy import deepcopy as dc
import random
fwd_table={'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}


backtable = {}
for key in fwd_table:
    val = fwd_table[key]
    if not val in backtable:
        backtable[val] = [key]
    else:
        backtable[val].append(key)
codons_not_to_mutate = []

class Triplet:
    def __init__(self,trip, verbose = False):
        self.verbose = verbose
        self.trip = trip.upper()
        assert(len(self.trip)==3) , 'Triplet not of length 3.'
        self.AA = fwd_table[self.trip]
        self.altcodons = dc(backtable[self.AA])
        if self.trip in self.altcodons and len(self.altcodons)>1: self.altcodons.remove(self.trip)
        self.orthoaltcodons = find_orthos(self.trip,self.altcodons)
        if self.verbose: print ('Triplet of %s initiated with amino acid being %s and alternative codons %s and %s orthogonal codons' %(self.trip,self.AA,self.altcodons,self.orthoaltcodons))

    def to_mutate(self):
        if self.trip in codons_not_to_mutate:
            return False
        else:
            return True

    def find_orthos(intrip,inlist):
        triscores =[]
        for tri in inlist:
            score = 0
            for x,c in enumerate(tri):
                if c == intrip[x]:
                    score += 1
            triscores.append(score)
        #print (triscores)
        minscore = min(triscores)
        outlist = []
        for x in range(0,len(inlist)):
            if triscores[x] == minscore:
                outlist.append(inlist[x])
        return outlist


    def random_alt(self):
        if self.to_mutate():
            random.seed()
            rannum = random.randint(0,len(self.altcodons)-1)
            return self.altcodons[rannum]
        else:
            return self.trip

def split_seq_to_Triplets(sequence, verbose=False):
    assert(len(sequence)%3 ==0), 'Sequence length not divisible by 3. Not gonna work with codon triplets.'
    list_of_triplet_objects = []
    for i in range(0,len(sequence),3):
        list_of_triplet_objects.append(Triplet(sequence[i:i+3],verbose))
    return list_of_triplet_objects

def new_random_sequence(sequence, verbose = False):
    triplets = split_seq_to_Triplets(sequence, verbose)
    newseq = ''.join([x.random_alt() for x in triplets])
    if verbose: print ('Writing old sequence %s to new sequence %s' %(sequence,newseq))
    return newseq

def simple_compare(seq1,seq2):
    assert(len(seq1)==len(seq2)), 'Sequence 1 not of same length as sequence 2. Obviously makes no sense like this.'
    hits = 0
    for i,c in enumerate(seq2):
        if c.lower() == seq1[i].lower():
            hits += 1
    percent = round(hits/len(seq1),2)
    return hits, percent
