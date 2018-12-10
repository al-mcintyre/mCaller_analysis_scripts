import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.stats import spearmanr,pearsonr

df = pd.read_csv('elisa_results.tab', sep='\s+')
colours = ['#857BB4','#E18D7D','#6ea5ff','#0099AD','#c2095a','#090446','#F71735','#f78c35','#951383','#60c9d1','#4c2faa','#feb95f']
#['#857BB4','#E18D7D','#65AAC3','#231F20','#BB4430','#773344','#79B791','#FF7E6B','#FFA69E','#DBD56E','#8D5A97','#FBBA72']
sns.set_palette(colours)
sns.set_style('white')

subdf = df.dropna(axis=0,how='any')
plt.figure(figsize=(4,4))
sns.regplot('ng_m6A','absorbance',data=subdf[subdf.date=='2017-12-14'],label='1')
sns.regplot('ng_m6A','absorbance',data=subdf[subdf.date=='2017-12-16'],label='2')
plt.xlabel('ng m6A')
plt.ylabel('OD')
plt.legend(title='replicate')
plt.title('Controls')
plt.tight_layout()
plt.savefig('elisa_controls.pdf',dpi=500,bbox_inches='tight')
plt.close()

expdf = []
for i,date in enumerate(['2017-12-14','2017-12-16']):
    pcdf = subdf[subdf.date==date]
    nc = pcdf[pcdf.species=='NC']['absorbance'].mean()
    pcdf["subtracted"] = pcdf.absorbance - nc
    pcdf = pcdf[pcdf.species != 'NC']
    #print pcdf.subtracted, nc
    slope, intercept, r_value, p_value, std_err = linregress(pcdf.ng_m6A,pcdf.subtracted)
    expdf.append(df[df.date==date])
    expdf[i] = expdf[i][pd.isnull(expdf[i]['ng_m6A'])] 
    expdf[i]["percent"] = (expdf[i].absorbance - nc)/2

xdf = pd.concat(expdf)
col_dict = {}
#print len(xdf[xdf.replicate == 1].percent),len(xdf[xdf.replicate == 2].percent)
#print xdf[xdf.replicate == 1].head(10), xdf[xdf.replicate == 2].head(10)
plt.figure(figsize=(4,4))
for i,species in enumerate(sorted(list(set(xdf["species"])))):
    sdf = xdf[xdf.species == species]
    plt.scatter(sdf[sdf.replicate == 1].percent,sdf[sdf.replicate == 2].percent,color=colours[i+2],label=species)
    col_dict[species] = colours[i+2]
fit = np.polyfit(xdf[xdf.replicate == 1].percent,xdf[xdf.replicate == 2].percent, deg=1)
plt.plot(xdf[xdf.replicate == 1].percent, fit[0] * np.asarray(xdf[xdf.replicate == 1].percent) + fit[1], color='#8997A1',alpha=0.3)
plt.xlabel('% m6A (rep 1)')
plt.ylabel('% m6A (rep 2)')
plt.legend()
plt.title('ELISA % m6A by species')
plt.tight_layout()
plt.savefig('elisa_replicate_comparison.pdf',dpi = 500,bbox_inches='tight')
plt.close()
rho,p = spearmanr(xdf[xdf.replicate == 1].percent,xdf[xdf.replicate == 2].percent)
print 'ELISA replicate comparison'
print rho, p
print pearsonr(xdf[xdf.replicate == 1].percent,xdf[xdf.replicate == 2].percent)

pb = pd.read_csv('pb_percent_m6A.tab', sep='\s+')
pbdf = pd.merge(xdf,pb,on='species')
plt.figure(figsize=(4,4))
pbx = []
elisay = []
for i,species in enumerate(sorted(list(set(pbdf["species"])))):
    sdf = pbdf[pbdf.species == species]
    plt.scatter(sdf[sdf.replicate == 1].pb_percent_m6Aover2xLEN,sdf["percent"].mean(),color=col_dict[species],label=species) #divide by 2 for old results
    pbx.append(sdf[sdf.replicate == 1].pb_percent_m6Aover2xLEN.mean())
    elisay.append(sdf["percent"].mean())
fit = np.polyfit(pbx, elisay, deg=1)
plt.plot(pbx, fit[0] * np.asarray(pbx) + fit[1], color='#8997A1',alpha=0.3)
plt.xlabel('% m6A (PacBio motifs)')
plt.ylabel('% m6A (ELISA mean)')
plt.legend()
plt.title('Motifs vs. ELISA')
plt.tight_layout()
plt.savefig('elisa_motif_comparison.pdf',dpi = 500,bbox_inches='tight')
plt.close()
rho,p = spearmanr(pbx,elisay)
print 'PB ELISA comparison'
print rho, p, 'Spearman rho and p'
print pearsonr(pbx,elisay)
