import sys
import re
import random
from Bio import SeqIO

species = 'ecoli'
fasta = '/athena/masonlab/scratch/projects/nasa/biomolecule_sequencer/data/pacbio/pb_ecoli_polished_assembly.fasta'
#grep m6A pb_ecoli_base_mods_4smrt_cells.gff > ecoli_k12_m6A_new.gff
#sed -i 's/scf7180000000002\|quiver/ecoli/g' ecoli_k12_m6A_new.gff 
#cut -d'|' -f 2- ecoli_k12_m6A_new.gff > ecoli_k12_m6A.gff 
gff='ecoli_k12_m6A.gff' #'/athena/masonlab/scratch/projects/nasa/biomolecule_sequencer/data/pacbio/pb_ecoli_m6A.gff' #'/athena/masonlab/scratch/projects/nasa/biomolecule_sequencer/data/pacbio/pb_ecoli_base_mods_4smrt_cells.gff'
#outfi = species+'_canonical_positions.txt'
training = species + '_training.txt'
testing = species + '_testing.txt'
unval_outfi = species + '_positions_to_comp_ecoli.txt'
pos_dict = {}

base_comps = {'A':'T','C':'G','T':'A','G':'C','N':'N','M':'M'}
motifs = {'GATC':'GATC','GCAC[ACGT]{6}GTT':'GCACNNNNNNGTT','AAC[ACGT]{6}GTGC':'AACNNNNNNGTGC'} 
motif_pos = {'GATC':1,'GCAC[ACGT]{6}GTT':2,'AAC[ACGT]{6}GTGC':1}

def comp(seq,base_comps=base_comps):
   return ''.join([base_comps[nt] for nt in list(seq)])

def revcomp(seq,strand='-'):
   if strand == '+':
      return seq
   else:
      return ''.join(list(comp(seq.upper()))[::-1])

strand2offset = {'-':-1,'+':1}

def check_motif(context,motifs=motifs,motif_pos=motif_pos):
    #print context, context[18:33], context[17:31]
    for motif in motifs:
        if re.match(motif,context[20-motif_pos[motif]:20+len(motifs[motif])]):
            return True
    return False

def find_motif(context,motifs=motifs,motif_pos=motif_pos):
    mos = []
    for motif in motifs:
        if re.match(motif,context[20-motif_pos[motif]:20+len(motifs[motif])]):
            mos.append(motifs[motif])
    if len(mos) > 0:
        return ','.join(mos)
    return 'none'

known_motifs,unknown_motifs = [],[]

incomplete_pos = set([])
pos_list = {}
lambda_inds = set(range(423317,471835)+range(698517,701053)) #based on BLASTn alignment
with open(gff,'r') as infi:
    for line in infi:
        if line[0] == '#':
            continue
        pos = int(line.split('\t')[3])-1 #adjust gff coordinates to 0-based
        strand = line.split('\t')[6]
        context = re.search('context=(.*?);',line).group(1)
        if context[20] != 'A':
            continue
        contig = line.split()[0]
        entry = (contig,pos,pos+1,strand)
        if (contig,strand) not in pos_list:
            pos_list[(contig,strand)] = []
        pos_list[(contig,strand)].append(pos)
        ##frac and QV score no provided in gff only if specified in SMRT analysis parameters
        ipd = re.search('IPDRatio=(.*?);',line).group(1)
        pos_dict[entry] = {'context':context,'IPD':ipd}
        try: 
            frac = float(re.search('frac=(.*?);',line).group(1))
            qv = float(line.split('identificationQv=')[1])
            if check_motif(context) and frac > 0.90 and qv >= 20 and pos not in lambda_inds:
                motif = find_motif(context)
                known_motifs.append(entry+(str(ipd),'m6A',context,'canon',motif))
            elif check_motif(context):
                motif = find_motif(context)
                unknown_motifs.append(entry+(str(ipd),'m6A',context,'canon_lowqual',motif)) 
            else:
                motif = 'none'
                unknown_motifs.append(entry+(str(ipd),'m6A',context,'noncanon',motif))
        except:
            motif = 'none'
            unknown_motifs.append(entry+(str(ipd),'m6A',context,'noncanon',motif))

controls = []
unmarked_canon = []
patterns = motifs.keys()
totals = {pat:[0,0] for pat in patterns}
found = 0
seqs= {}
for record in SeqIO.parse(fasta, "fasta"):
    contig = record.id
    seqs[contig] = str(record.seq).upper() 
    halfpoint = len(seqs[contig])/2

    #append all As in all contigs to potential controls
    controls = controls + [(contig,i,i+1,'+') for i, ltr in enumerate(seqs[contig]) if ltr == 'A' and i not in lambda_inds]
    controls = controls + [(contig,i,i+1,'-') for i, ltr in enumerate(seqs[contig]) if ltr == 'T' and i not in lambda_inds]

    #check whether canonical motifs in both strands found in pb gff of m6A positions
    for strand in ['+','-']:
        for pattern,offset in zip(patterns,[motif_pos[motif] for motif in patterns]):
            sq = revcomp(seqs[contig],strand)
            matches = re.finditer(pattern,sq)
            os = strand2offset[strand]*offset
            for match in matches:
                if strand == '-':
                    start = len(seqs[contig])-match.start(0)+os-1
                else:
                    start = match.start(0)+os
                if start in lambda_inds:
                    continue
                entry = (contig,start,start+1,strand)
                context = revcomp(seqs[contig][start-20:start+20],strand)
                totals[pattern][1] += 1
                if entry not in pos_dict and entry not in incomplete_pos: 
                    motif = motifs[pattern]
                    pos_dict[entry] = 1
                    #motif = find_motif(context)
                    unmarked_canon.append(entry+(0,'A',context,'missed',motif))
                else:
                    found +=1
                    totals[pattern][0] += 1 

print len(unmarked_canon),'canonical sites not found by PacBio,',found,'found'
for pat in totals:
    print motifs[pat]+':',totals[pat][0],'positions out of',totals[pat][1],'found:',totals[pat][0]*100./totals[pat][1],"%"

#make set of controls not included in pb m6A gff at canonical motif or noncanonical sites > 80%
def recursive_control_check(controls):
    choice = random.choice(controls)
    if choice in pos_dict or choice in incomplete_pos:
        return recursive_control_check(controls)
    else:
        pos_dict[choice] = 1
        return choice

ctrls = []
for x in range(max(len(known_motifs),len(unknown_motifs))):
    ctrl = recursive_control_check(controls)
    seq = seqs[ctrl[0]][ctrl[1]-20:ctrl[1]+21]
    if ctrl[3] == '-':
        seq = revcomp(seq)
    ctrls.append(ctrl+('0','A',seq,'control','none'))

#write bed files of positions
out = open(training,'w')
out2 = open(testing,'w')
print halfpoint
for pos in known_motifs:
    if int(pos[1]) <= halfpoint: 
        out.write('\t'.join([str(x) for x in pos])+'\n')
    else:
        out2.write('\t'.join([str(x) for x in pos])+'\n')
for pos in ctrls:
    if pos[1] <= halfpoint:
        out.write('\t'.join([str(x) for x in pos])+'\n')
    else:
        out2.write('\t'.join([str(x) for x in pos])+'\n')

print len(unknown_motifs),'noncanonical motifs'
print len(known_motifs),'known motifs' #,known_motifs[:10]
print len(unmarked_canon),'motifs missed'
print len(ctrls),'controls' #,ctrls[:10]
with open(unval_outfi,'w') as out:
    out.write('\t'.join(['species','position','pos+1','strand','IPD','base','context','type','motif'])+'\n')
    for pos in unknown_motifs:
        out.write('\t'.join([str(x) for x in pos])+'\n')
    for pos in known_motifs: #random.sample(known_motifs,min([len(unknown_motifs),len(known_motifs)])):
        out.write('\t'.join([str(x) for x in pos])+'\n')
    for pos in ctrls: #random.sample(ctrls,min([len(unknown_motifs),len(known_motifs)])):
        out.write('\t'.join([str(x) for x in pos])+'\n')
    for pos in unmarked_canon:
        out.write('\t'.join([str(x) for x in pos])+'\n')
