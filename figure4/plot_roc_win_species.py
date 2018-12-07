import sys 
import random
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

sns.set()
Gspecies = 'Ecoli'
extra = '/../20171212_8_zymo'
mcaller_fi = '..{}/{}/{}.methylation.positions.summary.bed'.format(extra,Gspecies,Gspecies.lower())
tombo_fi = '../{}_tombo/tombo_results.fraction_modified_reads'.format(Gspecies)
positions_fi = '..{}/{}/{}_positions_to_check.txt'.format(extra,Gspecies,Gspecies.lower())

positions = {'control':{},'canon':{}}

with open(positions_fi,'r') as infi1:
    for line in infi1:
        tab = line.strip().split('\t')
        chrom, pos, strand, sitetype, motif = tab[0],tab[2],tab[3], tab[7], tab[8]
        if sitetype in positions:
            positions[sitetype][(chrom,pos,strand)] = motif

mcaller_tp = {}
mcaller_tn = []
with open(mcaller_fi,'r') as infi2:
    for line in infi2:
        tab = line.strip().split('\t')
        #efaecalis  2421834 2421835 TCTTCMGGMMT 0.244791666667  +   192 5.074   37.407
        chrom,pos,strand,percent_meth = tab[0],tab[2],tab[5],float(tab[4])
        coords = (chrom,pos,strand) 
        if coords in positions['control']:
            mcaller_tn.append(percent_meth)
        elif coords in positions['canon']:
            motif = positions['canon'][coords]
            if motif not in mcaller_tp:
                mcaller_tp[motif] = []
            mcaller_tp[motif].append(percent_meth)

tombo_tp = {}
tombo_tn = []
strands = {'plus':'+','minus':'-'}
for s in strands:
    found = False
    with open('{}.{}.wig'.format(tombo_fi,s),'r') as infi3:
        for line in infi3:
            tab = line.strip().split()
            if tab[0] == 'variableStep':
                chrom = tab[1].split('=')[1]
                found = True
            elif found:
                pos,percent_meth = tab[0],float(tab[1])
                coords = (chrom,pos,strands[s]) 
                if coords in positions['control']:
                    tombo_tn.append(percent_meth)
                elif coords in positions['canon']:
                    motif = positions['canon'][coords]
                    if motif not in tombo_tp:
                        tombo_tp[motif] = []
                    tombo_tp[motif].append(percent_meth)

fig = plt.figure(figsize=(4,4))
#print mcaller_tp
print len(mcaller_tn), len(mcaller_tp)
pal1 = itertools.cycle(sns.color_palette("OrRd_d", len(tombo_tp)))
pal2 = itertools.cycle(sns.color_palette("PuBuGn_d", len(mcaller_tp)))
for motif in mcaller_tp:
    if len(motif.split(',')) == 1:
        tn = sorted(random.sample(mcaller_tn,min([len(mcaller_tp[motif]),len(mcaller_tn)])))
        tp = sorted(mcaller_tp[motif])
        preds = np.asarray( tp + tn )
        y_test = np.asarray( [1]*len(tp) + [0]*len(tn) , dtype=float)
        fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
        plt.plot(fpr,tpr,label='{} {}'.format('mCaller',motif), c=next(pal2))
for motif in tombo_tp:
    if len(motif.split(',')) == 1:
        tn = sorted(random.sample(tombo_tn,min([len(tombo_tp[motif]),len(tombo_tn)])))
        tp = sorted(tombo_tp[motif])
        preds = np.asarray( tp + tn )
        y_test = np.asarray( [1]*len(tp) + [0]*len(tn) , dtype=float)
        fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
        plt.plot(fpr,tpr,label='{} {}'.format('Tombo',motif),c=next(pal1))
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
plt.tight_layout()
plt.savefig('{}_roc_species.pdf'.format(Gspecies.lower()),bbox_inches='tight')
