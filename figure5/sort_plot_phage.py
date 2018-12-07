from Bio import SeqIO
import sys
import re
from scipy.stats import wilcoxon
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from markov_model import calc_markov_expected,calc_kr_expected,sort_phage,calc_ratio
from plotting_lib import *

iupac2regex = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
mtase_target_dict_iupac = { \
        'Enterococcus faecalis':['CRAANNNNNNRTTG','CTKVAG','CTCCAG'],\
        'Escherichia coli':['AAGANNNNNNNCTC','GATC'],\
        'Listeria monocytogenes':['GCANNNNNNNTGC'],\
        'Salmonella enterica':['CAGAG','BATGCATV','GATC'],\
        'Bacillus subtilis':['CNCANNNNNNNRTGT'],\
        'Staphylococcus aureus':['GAAGNNNNNNTTRG','TCTANNNNNNTTAA','GATCGVNY']\
        # no P aeruginosa sites
        }

colours = {'Bacillus subtilis':'#6ea5ff','Cryptococcus neoformans':'#0099AD','Escherichia coli':'#c2095a','Enterococcus faecalis':'#090446','Lactobacillus fermentum':'#F71735','Listeria monocytogenes':'#f78c35','Pseudomonas aeruginosa':'#951383','Staphylococcus aureus':'#60c9d1','Saccharomyces cerevisiae':'#4c2faa','Salmonella enterica':'#feb95f'}

if len(sys.argv[1:]) < 2:
    host = mtase_target_dict_iupac.keys()
    outfi = 'phage_all.fasta'
else:
    host = [' '.join(sys.argv[1:])]
    outfi = 'phage_' + host[0][0] + sys.argv[2] + '.fasta'

ratios = {}
seq_type = 'phage' #genome'
method = 'kr' #markov'
if seq_type == 'phage':
    virushost_fna = '../../REFS/viruses/virushostdb.genomic.fna'
    ids,control_ids,virus_types = sort_phage(host,mtase_target_dict_iupac)

comp = {'A':'T','C':'G','G':'C','T':'A','N':'N','K':'M','W':'S','M':'K','Y':'R','R':'Y','B':'A','D':'C','H':'G','V':'T','S':'W'}
for bacteria in host:
    Gspecies = bacteria[0]+bacteria.split()[1]
    seqs = []
    ratios[bacteria] = {}
    observed_hits = []
    shuffled_motif_hits = []
    markov_hits = []
    sequences,revcomps = [],[]
    if seq_type == 'prophage':
        virushost_fna = 'phaster/phaster_{}/phage_regions.fna'.format(Gspecies)
    if seq_type == 'genome':
        virushost_fna = '../{}/{}_pb.fasta'.format(Gspecies,Gspecies.lower())
    for record in SeqIO.parse(virushost_fna,"fasta"):
        if (seq_type == 'phage' and record.name in ids[bacteria]) or seq_type != 'phage':  
            seqs.append(record)
            sequence = str(record.seq).upper()
            sequences.append(sequence)
            revcomps.append(''.join([comp[N] for N in sequence[::-1]]))
            if seq_type != 'genome':
                markov_hits,observed_hits,ratios = calc_ratio(sequences,revcomps,bacteria,markov_hits,observed_hits,ratios,mtase_target_dict_iupac,method)
                #calc_ratio(sequences,revcomps,bacteria,markov_hits,observed_hits,ratios,mtase_target_dict_iupac,'markov')
                #break
    if seq_type == 'genome':
        markov_hits,observed_hits,ratios = calc_ratio(sequences,revcomps,bacteria,markov_hits,observed_hits,ratios,mtase_target_dict_iupac,method)

    print bacteria,len(seqs),'sequences found'
    print method,'motif enrichment:', sum(observed_hits),seq_type+' hits,',sum(markov_hits),method,'expected,',wilcoxon(observed_hits,markov_hits)

    #sys.exit(0)
    #plot results for species - moved to plotting_lib

    with open(outfi, "w") as output_handle:
        SeqIO.write(seqs, output_handle, "fasta")

all_ratios = [x for b in ratios for m in ratios[b] for x in ratios[b][m]]
#print all_ratios
plt.figure(figsize=(3,3))
sns.set_palette(['#7d1d3f'])
sns.distplot(all_ratios)
plt.plot([1,1],[0,1.3],linestyle='--')
plt.xlim([0,max(all_ratios)+1])
plt.ylim([0,1.3])
plt.xlabel('Ratio (observed/expected sites)')
plt.ylabel('Density')
plt.tight_layout()
plt.savefig('{}_motif_{}_representation.pdf'.format(seq_type,method),dpi=500,bbox_inches='tight')
plt.close()

if seq_type != 'genome':
    plot_phage(ratios,colours,seq_type,method)
else:
    plot_genome(ratios,colours,method)
