import sys
from copy import deepcopy as dc
import random

inputfile = sys.argv[1]
fin = open(inputfile,'r')
inputsequence = ''
for line in fin:
     inputsequence = inputsequence + line.strip('\n')
     inputsequence = inputsequence.replace(' ', '')

fwd_table={'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'TAA': '.', 'TAG': '.', 'TGA': '.'}

backtable = {}
for key in fwd_table:
    val = fwd_table[key]
    if not val in backtable:
        backtable[val] = [key]
    else:
        backtable[val].append(key)

exp_switch_backtable = dc(backtable)
exp_switch_backtable['R'].extend(backtable['K'])
exp_switch_backtable['K'].extend(backtable['R'])
exp_switch_backtable['D'].extend(backtable['Q'])
exp_switch_backtable['E'].extend(backtable['N'])
exp_switch_backtable['S'].extend(backtable['D'])
exp_switch_backtable['D'].extend(backtable['S'])

bad_codons = ['CTA', 'ATA', 'TAG', 'CGA', 'AGA', 'AGG', 'GGA']
surface = [15, 19, 20, 38, 41, 44, 45, 47, 50, 56, 59, 60, 63, 70, 79, 86, 91,106,109]
codons_not_to_mutate = []

class Triplet:
    def __init__(self,trip, verbose = False):
        self.verbose = verbose
        self.trip = trip.upper()
        assert(len(self.trip)==3) , 'Triplet not of length 3.'
        self.AA = fwd_table[self.trip]
        self.altcodons = dc(backtable[self.AA])
        for item in self.altcodons:
            if item in bad_codons: self.altcodons.remove(item)
        if self.trip in self.altcodons and len(self.altcodons)>1: self.altcodons.remove(self.trip)
        self.expaltcodons = dc(exp_switch_backtable[self.AA])
        for item in self.expaltcodons:
            if item in bad_codons: self.expaltcodons.remove(item)
        if self.trip in self.expaltcodons and len(self.expaltcodons)>1: self.expaltcodons.remove(self.trip)
        self.orthoaltcodons = self.find_orthos(self.trip,self.altcodons)
        self.orthoexpaltcodons = self.find_orthos(self.trip,self.expaltcodons)
        if self.verbose: print ('Triplet of %s initiated with amino acid being %s , alternative codons %s , alternative expanded codons %s , %s orthogonal codons and expanded orthogonal codons %s' %(self.trip,self.AA,self.altcodons,self.expaltcodons, self.orthoaltcodons, self.orthoexpaltcodons))

    def to_mutate(self):
        if self.trip in codons_not_to_mutate:
            return False
        else:
            return True

    def find_orthos(self,intrip,inlist):
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

    def random_expalt(self):
        if self.to_mutate():
            random.seed()
            rannum = random.randint(0,len(self.expaltcodons)-1)
            return self.expaltcodons[rannum]
        else:
            return self.trip

    def random_orthoalt(self):
        if self.to_mutate():
            random.seed()
            rannum = random.randint(0,len(self.orthoaltcodons)-1)
            return self.orthoaltcodons[rannum]
        else:
            return self.trip

    def random_orthoexpalt(self):
        if self.to_mutate():
            random.seed()
            rannum = random.randint(0,len(self.orthoexpaltcodons)-1)
            return self.orthoexpaltcodons[rannum]
        else:
            return self.trip

class Sequence:
    def __init__(self,inseq, verbose=False):
        self.inseq = inseq
        self.tripletlist = self.split_seq_to_Triplets()
        self.AAseq = ''
        for item in self.tripletlist:
            self.AAseq = self.AAseq + item.AA
        if verbose:
            print ('Sequence object initiated with input sequence:')
            print (self.inseq)
            print ('Amino acid sequence would be:')
            print (self.AAseq)

    def split_seq_to_Triplets(self, verbose=False):
        assert(len(self.inseq)%3 ==0), 'Sequence length not divisible by 3. Not gonna work with codon triplets.'
        list_of_triplet_objects = []
        for i in range(0,len(self.inseq),3):
            list_of_triplet_objects.append(Triplet(self.inseq[i:i+3],verbose))
        return list_of_triplet_objects

    def new_random_sequence(self, verbose = False):
        newseq = ''
        for i in range(0, len(self.tripletlist)):
            newseq = newseq + self.tripletlist[i].random_alt()
        if verbose: print ('Writing old sequence %s to new sequence %s' %(self.inseq,newseq))
        return newseq

    def new_random_exp_sequence(self, verbose = False):
        newseq = ''
        for i in range(0, len(self.tripletlist)):
            if i+2 in surface:
                newseq = newseq + self.tripletlist[i].random_expalt()
            else:
                newseq = newseq + self.tripletlist[i].random_alt()
        if verbose: print ('Writing old sequence %s to new sequence %s' %(self.inseq,newseq))
        return newseq

    def new_random_ortho_exp_sequence(self, verbose = False):
        newseq = ''
        for i in range(0, len(self.tripletlist)):
            if i+2 in surface:
                newseq = newseq + self.tripletlist[i].random_orthoexpalt()
            else:
                newseq = newseq + self.tripletlist[i].random_orthoalt()
        if verbose: print ('Writing old sequence %s to new sequence %s' %(self.inseq,newseq))
        return newseq

def simple_compare(seq1,seq2):
    assert(len(seq1)==len(seq2)), 'Sequence 1 not of same length as sequence 2. Obviously makes no sense like this.'
    hits = 0
    for i,c in enumerate(seq2):
        if c.lower() == seq1[i].lower():
            hits += 1
    percent = round(hits/len(seq1),2)
    return hits, percent

def graphical_compare(seq1,seq2):
    assert(len(seq1)==len(seq2)), 'Sequence 1 not of same length as sequence 2. Obviously makes no sense like this.'
    hits = 0
    for i,c in enumerate(seq2):
        if c.lower() == seq1[i].lower():
            hits += 1
    percent = round(hits/len(seq1),2)
    print ('Identical positions: %s. Perecent identical: %s.' %(hits,percent ))
    comparestring = ''
    for i in range(0,len(seq1)):
        if seq1[i].lower() == seq2[i].lower():
            comparestring = comparestring + '|'
        else:
            comparestring = comparestring + ' '
    end = 100
    print ('Showing only first %s positions.' %(end))
    print (seq1[:end])
    print (comparestring[:end])
    print (seq2[:end])
    return hits, percent
#newseq = ortho_sequences(inputsequence)

seqo = Sequence(inputsequence, True)

for i in range(0,10):
    newsequence = seqo.new_random_ortho_exp_sequence()
    newseqobj = Sequence(newsequence)
    newAA = newseqobj.AAseq
    Hits,Percent = simple_compare(seqo.inseq,newsequence)
    print ('New sequence with identical positions: %s. Percent identical: %s.' %(Hits, Percent))
    with open(inputfile.split('.')[-2]+'_new_'+str(i)+'.fasta','w') as fout:
        fout.write(seqo.new_random_ortho_exp_sequence())
    with open('new_AAsequence'+str(i)+'.fasta','w') as fout:
        fout.write('>new_AAsequence'+str(i)+'\n')
        fout.write(newAA)

'''Examples:
newsequence = seqo.new_random_ortho_exp_sequence()
graphical_compare(seqo.inseq,newsequence)



'''
