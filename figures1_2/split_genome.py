fi='/athena/masonlab/scratch/projects/nasa/biomolecule_sequencer/data/pacbio/pb_ecoli_polished_assembly.fasta'

fasta = open(fi,'r').read()
seq = ''.join(fasta.split('\n')[1:]).upper()
a_positions = [i for i, x in enumerate(list(seq)) if x == "A"]
t_positions = [i for i, x in enumerate(list(seq)) if x == "T"]
print len(seq), seq[:10]
with open('first_half.txt','w') as out:
    for pos in a_positions:
        print 
