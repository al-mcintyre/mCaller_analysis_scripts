base_comps = {'A':'T','C':'G','T':'A','G':'C','N':'N','M':'M'}

def comp(seq,base_comps=base_comps):
    return ''.join([base_comps[nt] for nt in list(seq)])

def revcomp(seq,strand='-'):
   if strand == '+':
      return seq.upper()
   else:
      return ''.join(list(comp(seq.upper()))[::-1])
