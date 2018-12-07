from Bio import SeqIO
import sys
import regex as re
from scipy.stats import wilcoxon
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from markov_model import calc_markov_expected,sort_phage,calc_ratio,calc_kr_expected
from plotting_lib import *

iupac2regex = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]','S':'[CG]','W':'[AT]','K':'[GT]','M':'[AC]'}
comp = {'A':'T','C':'G','G':'C','T':'A','N':'N','K':'M','W':'S','M':'K','Y':'R','R':'Y','B':'A','D':'C','H':'G','V':'T','S':'W'}

mtase_target_dict_iupac = {}
rm_type_dict = {}
typeset = set([])
motifs_per_rmtype = {}
infi = 'rusinov_RM_species_pairs.tsv'
method = 'markov'

host = set([])
for line in open(infi,'r'):
    if line[0] != '#':
        organism, motif, rm_type = line.split('\t')[2], line.split('\t')[3], line.split('\t')[6]
        species = ' '.join(organism.split()[:2])
        host.add(species)
        if species not in mtase_target_dict_iupac:
            mtase_target_dict_iupac[species] = []
        if motif not in set(mtase_target_dict_iupac[species]):
            mtase_target_dict_iupac[species].append(motif)
        if rm_type not in motifs_per_rmtype:
            motifs_per_rmtype[rm_type] = set([])
        rm_type_dict[(species,motif)] = rm_type
        motifs_per_rmtype[rm_type].add((species,motif))
        typeset.add(rm_type)

for rm_type in motifs_per_rmtype:
    print '{}: {} motifs'.format(rm_type,len(motifs_per_rmtype[rm_type]))
    if rm_type == 'Type IV':
        print motifs_per_rmtype[rm_type]

ratios = {}
virushost_fna = '../../REFS/viruses/virushostdb.genomic.fna'
print '{} host species found'.format(len(host))
ids,control_ids,virus_types = sort_phage(list(host),mtase_target_dict_iupac)

ratio_dict = {}
ratios = {}
pairs = {'host':[],'virus':[],'motif':[],'type':[],'max_markov_ratio':[]}
examples = 0
for bacteria in host:
    ratios[bacteria] = {}
    for record in SeqIO.parse(virushost_fna,"fasta"):
        if record.name in ids[bacteria]:
            sequences = [str(record.seq).upper()]
            revcomps = [''.join([comp[N] for N in sequences[0][::-1]])]
            for iupac_target in mtase_target_dict_iupac[bacteria]:
                mtase_target = ''.join([iupac2regex[i] for i in list(iupac_target)])
                Nx,Ex = 0,0
                for sequence,revcomp in zip(sequences,revcomps):
                    if method == 'markov' or method == 'markov1':
                        Nx += len(re.findall(mtase_target,sequence,overlapped=True)) + len(re.findall(mtase_target,revcomp,overlapped=True))
                        Ex += calc_markov_expected(iupac_target,sequence,True,method)
                        if Ex > 0:
                            ratio = Nx*1./Ex
                        else:
                            ratio = (Nx+1)*1./(Ex+1)
                    elif method == 'kr':
                        ratio = calc_kr_expected(iupac_target,sequence,rev=False)
                rm_type = rm_type_dict[(bacteria,iupac_target)] 
                if (method != 'kr') or (method == 'kr'): # and Ex > 0) or (method == 'kr'):
                    if rm_type not in ratio_dict:
                           ratio_dict[rm_type] = {}
                    if iupac_target not in ratio_dict[rm_type]:
                        ratio_dict[rm_type][iupac_target] = []
                        examples +=1
                    ratio_dict[rm_type][iupac_target].append(ratio)
                    pairs['host'].append(bacteria)
                    pairs['virus'].append(record.name)
                    pairs['motif'].append(iupac_target)
                    pairs['type'].append(rm_type)
                    pairs['max_markov_ratio'].append(ratio)
                    if examples%100 == 0:
                        print examples
                        #break
        #if examples == 10:
         #   break
    #if examples == 10:
    #    break

for rm_type in ratio_dict:
    print '{}: {} phage genomes'.format(rm_type,len(ratio_dict[rm_type]))

pairs_df = pd.DataFrame(pairs)
pairs_df = pairs_df[['type','motif','host','virus','max_markov_ratio']]
pairs_df = pairs_df.sort_values(by=['type','motif','host','virus'])
pairs_df.to_csv('motif_phage_table.tsv',sep='\t',index=False, float_format='%.3f')
plot_density_curves(ratio_dict,method)
