import sys
import matplotlib as ml
from matplotlib import pyplot as pp
import numpy as np


if len(sys.argv)>0:
    infile = sys.argv[1]
else:
    print ('Please specify inputfile with sequence in commandline.')
    exit()

windowsize = 30
if len(sys.argv)>1:
    windowsize = sys.argv[2]


inputsequence = ''
with open(infile,'r') as fin:
    for line in fin:
         inputsequence = inputsequence + line.strip('\n')
         inputsequence = inputsequence.replace(' ', '')

positions = []
gccontent_list = []

for x in range(0,len(inputsequence)-windowsize):
    counter = 0
    for a in inputsequence[x:x+windowsize]:
        if a.lower() == 'g' or a.lower() == 'c':
            counter += 1
    percent = round(counter/windowsize*100,2))
    gccontent_list.append(percent)
    positions.append(x)

pos = np.array(positions)
gc = np.array(gccontent_list)
pp.plot(pos,gc)
pp.xlabel('Position')
pp.ylabel('GC content in %')
pp.title('GC Content with %s window of %s' %(windowsize,infile))
pp.show()
