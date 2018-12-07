import sys 
import random
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

sns.set()
Gspecies = 'Efaecalis' #'Ecoli' #faecalis'
extra = '' #/../20171212_8_zymo'
mcaller_fi = '..{}/{}/{}.methylation.positions.summary.bed'.format(extra,Gspecies,Gspecies.lower())
#tombo_results.fraction_m6A_reads.minus.wig
tombo_fi = '../{}_tombo/tombo_results.fraction_m6A_reads'.format(Gspecies)
positions_fi = '../{}/{}_positions_to_check.txt'.format(Gspecies,Gspecies.lower())
neg_control_fi = '../Paeruginosa/paeruginosa.methylation.{}.positions.summary.bed'.format(Gspecies.lower())
neg_control_positions = '../Paeruginosa/paeruginosa_positions_to_comp_{}.txt'.format(Gspecies.lower())
tombo_pseudo = '../Paeruginosa_tombo/tombo_results.fraction_modified_reads'

positions = {'control':{},'canon':{}}

with open(positions_fi,'r') as infi1:
    for line in infi1:
        tab = line.strip().split('\t')
        chrom, pos, strand, sitetype, motif = tab[0],tab[2],tab[3], tab[7], tab[8]
        if sitetype == 'canon': # in positions:
            positions[sitetype][(chrom,pos,strand)] = motif

with open(neg_control_positions,'r') as infi1:
    for line in infi1:
        tab = line.strip().split('\t')
        chrom, pos, strand, sitetype, motif = tab[0],tab[2],tab[3], tab[7], tab[8]
        if motif != 'none': # in positions:
            for mo in motif.split(','):
                positions['control'][(chrom,pos,strand)] = mo

mcaller_tp = {}
mcaller_tn = {}
with open(mcaller_fi,'r') as infi2:
    for line in infi2:
        tab = line.strip().split('\t')
        #efaecalis  2421834 2421835 TCTTCMGGMMT 0.244791666667  +   192 5.074   37.407
        chrom,pos,strand,percent_meth = tab[0],tab[2],tab[5],float(tab[4])
        coords = (chrom,pos,strand) 
        #if coords in positions['control']:
        #    mcaller_tn.append(percent_meth)
        if coords in positions['canon']:
            motif = positions['canon'][coords]
            if motif not in mcaller_tp:
                mcaller_tp[motif] = []
            mcaller_tp[motif].append(percent_meth)

for motif in mcaller_tp:
    print 'mCaller:', motif, len([p for p in mcaller_tp[motif] if p >= 0.5]), 'of', len(mcaller_tp[motif])

with open(neg_control_fi,'r') as infi2:
    for line in infi2:
        tab = line.strip().split('\t')
        chrom,pos,strand,percent_meth = tab[0],tab[2],tab[5],float(tab[4])
        coords = (chrom,pos,strand)
        if coords in positions['control']:
            motif = positions['control'][coords]
            if motif not in mcaller_tn:
                mcaller_tn[motif] = []
            mcaller_tn[motif].append(percent_meth)

tombo_tp = {}
tombo_tn = {}
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
                if coords in positions['canon']:
                    motif = positions['canon'][coords]
                    if motif not in tombo_tp:
                        tombo_tp[motif] = []
                    tombo_tp[motif].append(percent_meth)
    found = False
    with open('{}.{}.wig'.format(tombo_pseudo,s),'r') as infi3:
        for line in infi3:
            tab = line.strip().split()
            if tab[0] == 'variableStep':
                chrom = tab[1].split('=')[1]
                found = True
            elif found:
                pos,percent_meth = tab[0],float(tab[1])
                coords = (chrom,pos,strands[s])
                if coords in positions['control']:
                    motif = positions['control'][coords]
                    if motif not in tombo_tn:
                        tombo_tn[motif] = []
                    tombo_tn[motif].append(percent_meth)

for motif in tombo_tp:
    print 'Tombo:', motif, len([p for p in tombo_tp[motif] if p >= 0.5]), 'of', len(tombo_tp[motif])

fig = plt.figure(figsize=(7,4))
#print mcaller_tp
print len(mcaller_tn)
pal1 = itertools.cycle(sns.color_palette("OrRd_d", len(tombo_tp)))
pal2 = itertools.cycle(sns.color_palette("PuBuGn_d", len(mcaller_tp)))
for motif in mcaller_tp:
    if motif in mcaller_tn:
        tn = sorted(mcaller_tn[motif]) #random.sample(mcaller_tn,min([len(mcaller_tp[motif]),len(mcaller_tn)])))
        tp = sorted(mcaller_tp[motif])
        preds = np.asarray( tp + tn )
        y_test = np.asarray( [1]*len(tp) + [0]*len(tn) , dtype=float)
        print 'mCaller: {} true negatives, {} true positives for {}. AUC = {}'.format(len(tn),len(tp),motif,metrics.roc_auc_score(y_test,preds))
        min_true = min([len(tn),len(tp)])
        print len([1 for x in random.sample(tp,min_true) if x >= 0.5]), len([0 for x in random.sample(tn,min_true) if x < 0.5]), min_true, len(tn), len(tp)
        balanced_preds = [1 if x >= 0.5 else 0 for x in random.sample(tp,min_true)] + [0 if x < 0.5 else 1 for x in random.sample(tn,min_true)]
        print 'Accuracy = {}'.format(metrics.accuracy_score([1]*min_true + [0]*min_true,balanced_preds))
        fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
        plt.plot(fpr,tpr,label='{} {}'.format('mCaller',motif), c=next(pal2))
for motif in tombo_tp:
    if motif in mcaller_tn:
        tn = sorted(mcaller_tn[motif]) #random.sample(tombo_tn,min([len(tombo_tp[motif]),len(tombo_tn)])))
        tp = sorted(tombo_tp[motif])
        preds = np.asarray( tp + tn )
        y_test = np.asarray( [1]*len(tp) + [0]*len(tn) , dtype=float)
        fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
        print 'Tombo: {} true negatives, {} true positives for {}. AUC = {}'.format(len(tn),len(tp),motif,metrics.roc_auc_score(y_test,preds))
        plt.plot(fpr,tpr,label='{} {}'.format('Tombo',motif),c=next(pal1))
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
plt.tight_layout()
plt.savefig('{}_roc_motifs.pdf'.format(Gspecies.lower()),bbox_inches='tight')
