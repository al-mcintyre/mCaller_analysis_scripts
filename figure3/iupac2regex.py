def i2r(iupac,rev=False):
    iupac2regex = iupac2regex = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]','K':'[GT]','M':'[AC]','W':'[AT]','S':'[GC]'}
    if rev:
        comp = {'A':'T','C':'G','G':'C','T':'A','N':'N','K':'M','W':'S','M':'K','Y':'R','R':'Y','B':'A','D':'C','H':'G','V':'T','S':'W'}
        iupac_target = ''.join([comp[i] for i in list(iupac)])
    else:
        iupac_target = iupac
    return(''.join([iupac2regex[i] for i in list(iupac_target)]))
