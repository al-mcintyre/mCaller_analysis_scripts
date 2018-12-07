import sys
import re
import random
from Bio import SeqIO
from iupac2regex import i2r
from basic_seq import revcomp

class species_info:
    def __init__(self):
        self.Gdot = ''
        self.Gspecies = ''
        self.species = ''
        self.fasta = ''
        self.gff = ''
        self.outfi = ''
        self.motifs = {}
    def update(self,line):
        Gdot = line.split()[0]
        if self.Gdot != Gdot and Gdot != 'species':
            self.Gdot = Gdot
            self.Gspecies = ''.join(Gdot.split('.'))
            self.species = self.Gspecies.lower()
            self.fasta = '../{}/{}_pb.fasta'.format(self.Gspecies,self.species)
            self.gff = '../{}/{}_mm.gff'.format(self.Gspecies,self.species)
            self.outfi = '../{}/{}_positions_to_check.txt'.format(self.Gspecies,self.species)
            self.motifs = {}
    def add_motif(self,line):
        Gdot = line.split()[0]
        if Gdot != 'species':
            motif = line.split()[1]
            if len(motif.split('.')) > 1:
                motif = motif.split('.')[1]
            if motif != 'none':
                self.motifs[i2r(motif)] = {'iupac':motif,'pos':int(line.split()[2])}
    def new_species(self,line):
        Gdot = line.split()[0]
        if self.Gdot != Gdot and Gdot != 'species' and self.Gdot != '':
            return True
        else:
            return False

def find_motif(context,motifs):
    mos = []
    mlist = sorted(motifs.keys())
    for motif in mlist:
        if re.match(motif,context[20-motifs[motif]['pos']:20+len(motifs[motif]['iupac'])]):
            mos.append(motifs[motif]['iupac'])
    if len(mos) > 0:
        return ','.join(mos)
    return 'none'

def pacbio_results(gff,motifs,pos_list,pos_dict):
    pb_found,nonmotif_sites = [],[]
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
            entry = (contig,str(pos),str(pos+1),strand)
            if (contig,strand) not in pos_list:
                pos_list[(contig,strand)] = []
            pos_list[(contig,strand)].append(pos)
            ##frac and QV score no longer provided in gff by default
            ipd = re.search('IPDRatio=(.*?)\n',line).group(1)
            pos_dict[entry] = {'context':context,'IPD':ipd}
            motif = find_motif(context,motifs)
            if motif != 'none':
                pb_found.append(entry+(str(ipd),'m6A',context,'canon',motif))
            else: 
                nonmotif_sites.append(entry+(str(ipd),'m6A',context,'noncanon',motif))
    return pb_found,nonmotif_sites,pos_dict,pos_list

#make set of controls not included in pb m6A gff 
def recursive_control_check(controls,pos_dict):
    choice = random.choice(controls)
    if choice in pos_dict:
        return recursive_control_check(controls)
    else: 
        return choice

def get_controls(fasta,motifs,pos_dict,pb_found,nonmotif_sites,strand2offset):
    controls,pb_missed = [],[]
    found = 0
    seqs = {}
    totals = {motif:[0,0] for motif in motifs}
    for record in SeqIO.parse(fasta, "fasta"):
        contig = record.id
        seqs[contig] = str(record.seq).upper()

        #append all As in all contigs to potential controls
        controls = controls + [(contig,i,i+1,'+') for i, ltr in enumerate(seqs[contig]) if ltr == 'A']
        controls = controls + [(contig,i,i+1,'-') for i, ltr in enumerate(seqs[contig]) if ltr == 'T']
        #check whether canonical motifs in both strands found in pb gff of m6A positions
        for strand in ['+','-']:
            for motif in motifs:
                sq = revcomp(seqs[contig],strand)
                matches = re.finditer(motif,sq)
                os = strand2offset[strand]*motifs[motif]['pos']
                for match in matches:
                    if strand == '-':
                        start = len(seqs[contig])-match.start(0)+os-1
                    else:
                        start = match.start(0)+os
                    entry = (contig,str(start),str(start+1),strand)
                    if strand == '-':
                        start = start + 1
                    context = revcomp(seqs[contig][start-20:start+20],strand)
                    motifsInCtxt = find_motif(context,motifs)
                    totals[motif][1] += 1
                    if entry not in pos_dict:
                        full_entry = entry+(0,'A',context,'missed',motifsInCtxt)
                        if full_entry not in set(pb_missed):
                            pb_missed.append(full_entry)
                    else:
                        totals[motif][0] += 1
                        found += 1
    print len(pb_missed),'canonical sites not found by PacBio,',found,'found'
    for motif in totals:
        print '{}: {} positions out of {} found = {}'.format(motifs[motif]['iupac'],totals[motif][0],totals[motif][1],totals[motif][0]*1./totals[motif][1])

    ctrls = []
    for x in range(max(len(pb_found),len(nonmotif_sites))):
        control = recursive_control_check(controls,pos_dict)
        seq = seqs[control[0]][control[1]-20:control[1]+20]
        if control[2] == '-':
            seq = revcomp(seq)
        ctrls.append(control+('0','A',seq,'control','none'))
    return pb_missed,ctrls

def write_bed(outfi,pb_found,nonmotif_sites,pb_missed,controls):
    print len(nonmotif_sites),'noncanonical motifs'
    print len(pb_found),'known motifs'
    print len(pb_missed),'motifs missed'
    print len(controls),'controls'
    with open(outfi,'w') as out:
        out.write('\t'.join(['species','position','pos+1','strand','IPD','base','context','type','motif'])+'\n')
        for pos in nonmotif_sites:
            out.write('\t'.join([str(x) for x in pos])+'\n')
        for pos in pb_found: 
            out.write('\t'.join([str(x) for x in pos])+'\n')
        for pos in controls: 
            out.write('\t'.join([str(x) for x in pos])+'\n')
        for pos in pb_missed:
            out.write('\t'.join([str(x) for x in pos])+'\n')

def analyze_last_species(sinfo):
    pos_dict,pos_list = {},{}
    strand2offset = {'-':-1,'+':1}
    pb_found,nonmotif_sites,pos_dict,pos_list = pacbio_results(sinfo.gff,sinfo.motifs,pos_dict,pos_list)
    pb_missed,controls = get_controls(sinfo.fasta,sinfo.motifs,pos_dict,pb_found,nonmotif_sites,strand2offset)
    write_bed(sinfo.outfi,pb_found,nonmotif_sites,pb_missed,controls)

with open('motif_summary.txt','r') as motif_fi:
    sinfo = species_info()
    for line in motif_fi:
        if sinfo.new_species(line):
            print '\n'+sinfo.Gdot
            analyze_last_species(sinfo)
            #break
        sinfo.update(line)
        sinfo.add_motif(line)
    print '\n'+sinfo.Gdot
    analyze_last_species(sinfo)
