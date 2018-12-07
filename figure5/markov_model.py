from Bio import SeqIO
import sys
import regex as re
from scipy.stats import wilcoxon
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from itertools import combinations

iupac2regex = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]','S':'[CG]','W':'[AT]','K':'[GT]','M':'[AC]'}
comp = {'A':'T','C':'G','G':'C','T':'A','N':'N','K':'M','W':'S','M':'K','Y':'R','R':'Y','B':'A','D':'C','H':'G','V':'T','S':'W'} 

def calc_markov_expected(target,sequence,rev=True,order='markov'):
    ts = list(target)
    if order =='markov':
        ms = [len(target)-1,len(target)-2] #for maximal order, -1 and -2 
    else:
        ms = [order+1,order]
    numdenom = [[],[]]
    subsequences = [[],[]]
    if rev:
        try:
            revcomp = ''.join([comp[N] for N in sequence[::-1]])
        except KeyError:
            print sequence
            sys.exit(0)
    for n in range(2):
        #print len(ts)-ms[n]
        for i in range(len(ts)-ms[n]+1-2*n):
            targ = ''.join([iupac2regex[ts[j]] for j in range(n+i,n+i+ms[n])])
            if n == 1:
                targ = targ+'[ACGT]'
            subsequences[n].append(targ)
            #print n, targ, len(re.findall(targ,sequence))
            count = len(re.findall(targ,sequence,overlapped=True))
            if rev:
                count+=len(re.findall(targ,revcomp,overlapped=True))
            numdenom[n].append(count) 
    numcp = numdenom[0][:]
    for n in numcp:
        if n in set(numdenom[1]):
            numdenom[0].remove(n)
            numdenom[1].remove(n)
    divided = np.divide(numdenom[0][:-1],numdenom[1],dtype=float)
    prod = 1
    for d in divided:
        prod = prod*d
    expected = prod*numdenom[0][-1] #*(len(sequence)-len(target)+1)
    #var = expected*
    #print target+'\t',str(expected) #+'\n\t',subsequences[0],'\n\t',subsequences[1]
    #print '\t',numdenom[0],'\n\t',numdenom[1]
    return expected

def get_submotifs(target,length):
    submotifs = []
    list_of_positions = range(len(target))
    combos = combinations(list_of_positions, length)
    for combo in list(combos):
        new_motif = ''
        for i in list_of_positions:
            if i in set(combo):
                new_motif = new_motif+iupac2regex[target[i]]
            else:
                new_motif = new_motif+'[ACGT]'
        submotifs.append(new_motif)
    return(submotifs)

def calc_kr_expected(target,sequence,rev=True):
    ts = list(target)
    numdenom = [[],[]]
    subsequences = [[],[]]
    if rev:
        try:
            revcomp = ''.join([comp[N] for N in sequence[::-1]])
        except KeyError:
            print sequence
            sys.exit(0)
    for n in range(1,len(target)+1):
        lw_ls = (len(target)-n)%2 #ie. -1^(L(w)-L(s))
        submotifs = get_submotifs(target,n)
        for submotif in submotifs:
            count = len(re.findall(submotif,sequence,overlapped=True))
            if n == len(target):
                print submotif,count,np.divide(count,len(sequence),dtype=float)
            if rev:
                count+=len(re.findall(submotif,revcomp,overlapped=True))
            subsequences[lw_ls].append(submotif)
            numdenom[lw_ls].append(count) #np.divide(count,len(sequence),dtype=float))
    
    numcp = numdenom[0][:]
    for n in numcp:
        if n in set(numdenom[1]):
            numdenom[0].remove(n)
            numdenom[1].remove(n)
    minsubseq,maxsubseq = min(len(numdenom[0]),len(numdenom[1])), max(len(numdenom[0]),len(numdenom[1]))
    if len(numdenom[0]) < maxsubseq:
        numdenom[0] = numdenom[0]+[1]*(maxsubseq-minsubseq)
    else:
        numdenom[1] = numdenom[1]+[1]*(maxsubseq-minsubseq)
    #print minsubseq, len(numdenom[0]),len(numdenom[1])
    #print numdenom[0], numdenom[1]
    divided = np.divide(numdenom[0],numdenom[1],dtype=float)
    prod = 1
    if len(target)%2 == 0:
        prod = prod*len(sequence)
    else:
        prod = np.divide(prod,len(sequence),dtype=float)
    high_divided = sorted([x for x in divided if x >= 1])[::-1]
    low_divided = sorted([x for x in divided if x < 1])
    mindiv = min(len(high_divided),len(low_divided))
    for hd,ld in zip(high_divided[:mindiv],low_divided[:mindiv]):
        prod = prod*hd*ld
    for extra in high_divided[mindiv:]+low_divided[mindiv:]:
        prod = prod*extra
    kr = prod
    #print target+'\t',str(expected) #+'\n\t',subsequences[0],'\n\t',subsequences[1]
    #print ' '.join(subsequences[0]),'\n'
    #print ' '.join(subsequences[1]),'\n'
    #print '\t',numdenom[0],'\n\t',numdenom[1]
    return kr

def calc_ratio(sequences,revcomps,bacteria,markov_hits,observed_hits,ratios,mtase_target_dict_iupac,method):
    Nx,Ex = 0,0
    target_types = []
    for iupac_target in mtase_target_dict_iupac[bacteria]:
        if iupac_target not in ratios[bacteria]:
            ratios[bacteria][iupac_target] = []
        mtase_target = ''.join([iupac2regex[i] for i in list(iupac_target)])
        if method == 'markov':
            for sequence,revcomp in zip(sequences,revcomps):
                Nx += len(re.findall(mtase_target,sequence,overlapped=True)) + len(re.findall(mtase_target,revcomp,overlapped=True))
                Ex += calc_markov_expected(iupac_target,sequence)
            observed_hits.append(Nx)
            markov_hits.append(Ex)
            #Ex,Varx = calc_markov_expected(iupac_target,sequence)
            #Zx = (Nx-Ex)*1./Varx
            if bacteria == 'Bacillus subtilis':
                print 'markov = {}, observed = {}'.format(Ex, Nx)
            if markov_hits[-1] > 0:
                ratio = Nx*1./Ex
            else:
                ratio = (Nx+1)*1./(Ex+1)
            ratios[bacteria][iupac_target].append(ratio)
            if bacteria == 'Bacillus subtilis':
                print 'markov = {}, observed = {}, ratio = {}'.format(Ex, Nx, ratio)
        elif method == 'kr':
            ratio = calc_kr_expected(iupac_target,'NNNNNNNNNN'.join(sequences),rev=False)
            ratios[bacteria][iupac_target].append(ratio)
            print 'Kr = {}'.format(ratio)
    return markov_hits,observed_hits,ratios

def sort_phage(host,mtase_target_dict_iupac):
    virushost_tsv = '../../REFS/viruses/virushostdb.tsv'
    ids = {bacteria:set([]) for bacteria in host}
    control_ids = set([])
    virus_types = {}
    phage_count,control_count = 0,0
    for line in open(virushost_tsv,'r'):
        if len(line) > 1:
            virus_ids, virus_host_domain, virus_host = line.split('\t')[3], line.split('\t')[9].split(';')[0], ' ' .join(line.split('\t')[8].split()[:2])
            try:
                virus_type = line.split(';')[1].split(',')[0].strip()
                if (virus_host == host[0] or (len(host) > 1 and virus_host in mtase_target_dict_iupac)) and virus_type == 'dsDNA viruses': #sort phage by type
                    for virus_id in virus_ids.split(', '):
                        ids[virus_host].add(virus_id)
                        virus_types[virus_id] = virus_type
                        phage_count += 1
                if (virus_host_domain == 'Eukaryota'):
                    for virus_id in virus_ids.split(', '):
                        control_ids.add(virus_id)
                        virus_types[virus_id] = virus_type
                        control_count += 1
            except (IndexError,TypeError) as e:
                print line

    print '{} records found'.format(phage_count)
    print '{} control sequences found'.format(control_count)
    #for h in host:
    #    print '\t{} : {} phage'.format(h,len(ids[h]))
    return ids, control_ids, virus_types
