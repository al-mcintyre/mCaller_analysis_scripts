import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools

sns.set_style('whitegrid')
sns.set_palette(["#0b3c49","#5f9fb9"])

motif_summary = pd.read_csv('motif_summary.txt',sep='\t')
motif_summary['method'] = 'PB'
motif_summary2 = pd.read_csv('motif_summary.txt',sep='\t')
motif_summary2['method'] = 'ONT'
motif_summary3 = pd.read_csv('motif_summary.txt',sep='\t')
motif_summary3['method'] = 'MeDIP'
motif_summary3['fraction_m6A'] = motif_summary2['fraction_medip']
motif_summary2['fraction_m6A'] = motif_summary2['fraction_ONT']
motif_summary['fraction_m6A'] = motif_summary['fraction_PB']
motifs = pd.concat([motif_summary,motif_summary2,motif_summary3])

nonpb_motifs = pd.concat([motif_summary2,motif_summary3])
motifs_none = nonpb_motifs[nonpb_motifs.motif == 'none']
motifs_none['species_method'] = motifs_none['species'] + " " + motifs_none['method']
motifs_none.sort_values(by='species',inplace=True)

motifs = motifs[motifs.motif != 'none']

for df,name in zip([motifs,motifs_none],['motifs','nonmotifs']):
   if name == 'nonmotifs':
      fig = plt.figure(figsize=(4,4))
      #sns.set_palette(['#65AAC3','#BB4430','#773344','#FF7E6B','#FFA69E','#DBD56E','#8D5A97','#FBBA72'])
      sns.set(font_scale=1.3,style='whitegrid') 
      sns.set_palette(['#6ea5ff','#95bdff','#c2095a','#d24c87','#090446','#4c4878','#f78c35','#f9ab6c','#951383','#b153a4','#60c9d1','#8bd7dd','#4c2faa','#7c67c1','#feb95f','#fecc8a'])
      g = sns.catplot(x='motif',y='fraction_m6A',hue='species_method',data=df,kind='bar',aspect=1.4,legend=False,facet_kws={'legend_out':False})
      hatches = itertools.cycle(['','//'])
      for i, bar in enumerate(g.axes[0,0].patches):
         hatch = next(hatches)
         bar.set_hatch(hatch)
      g.axes[0,0].set_ylim([0,100])
      g.axes[0,0].set_xticklabels([''])
      g.axes[0,0].set_xlabel('nonmotif sites (PacBio)')
      g.axes[0,0].set_ylabel('percent of sites m6A\n(PacBio + mCaller or MeDIP-seq)')
   else:
      fig = plt.figure(figsize=(10,4))
      sns.set(font_scale=1.3,style='whitegrid')
      sns.set_palette(["#133741","#6a9aae","#a5c7c4"])
      g = sns.catplot(x='motif',y='fraction_m6A',hue='method',data=df,kind='bar',aspect=2.5,legend=False,facet_kws={'legend_out':False})
      g.set_xticklabels(rotation=90, fontsize=13)
      for ax in g.axes.flat:
         box = ax.get_position()
         ax.set_xlabel('')
         ax.set_position([box.x0,box.y0,box.width*0.85,box.height])
      plt.legend(loc='upper left',bbox_to_anchor=(1,0.5))
      g.axes[0,0].set_ylabel('percent of sites m6A')

   plt.tight_layout()

   plt.savefig(name+'_plot.pdf',dpi=800,bbox_inches='tight',transparent=True)
