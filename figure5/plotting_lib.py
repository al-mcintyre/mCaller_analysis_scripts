import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plot_phage(ratios,colours,seq_type,method):
    #plot by motif, coloured by species
    sns.set_palette(['#D3D3D3'])
    fig = plt.figure(1,figsize=(6,3))
    ax = fig.add_subplot(111)
    i = 1
    ylabs = []
    for bacteria in sorted(ratios.keys()):
        for motif in ratios[bacteria]:
            print motif, i
            ylabs.append(motif)
            bp = ax.boxplot(ratios[bacteria][motif],positions=range(i,i+1),notch=True,vert=False,patch_artist=True,widths=0.4)
            for box in bp['boxes']:
                box.set(color=colours[bacteria],linewidth=2,facecolor=colours[bacteria])
            plt.text(min(max(ratios[bacteria][motif])+0.2,2.5),i-0.2,'N = {}'.format(len(ratios[bacteria][motif])),fontsize=9)
            i+=1
        #sns.distplot(ratios[bacteria],hist=False,color=colours[bacteria])
    max_x = 3 #max(max(all_ratios)+max(all_ratios)*0.5,2)
    ax.set_yticks(range(1,len(ylabs)+1))
    ax.set_yticklabels(ylabs)
    plt.plot([1,1],[0,(len(ylabs)+1)],linestyle='--')
    plt.xlim([0,max_x])
    plt.ylim([0,(len(ylabs)+1)])
    plt.xlabel('Ratio (observed/expected sites)')
    plt.tight_layout()
    plt.savefig('{}_motif_box_by_species_{}.pdf'.format(seq_type,method),dpi=500,bbox_inches='tight',transparent=True)
    plt.close()

def plot_genome(ratios,colours,method,seq_type='genome'):
    sns.set_palette(['#D3D3D3'])
    fig = plt.figure(1,figsize=(6,1.4))
    ax = fig.add_subplot(111)
    for bacteria in ratios:
        for motif in ratios[bacteria]:
            ratio = ratios[bacteria][motif][0]
            print motif, ratio
            ax.plot([ratio,ratio],[-0.3,0.3],color=colours[bacteria])
            if len(ratios[bacteria]) > 1:
                plt.text(ratio-0.005,0.6,motif[0],fontsize=9) #label with first letter if multiple
    max_x = 1.2
    ax.set_yticks([])
    plt.plot([1,1],[-1,1],linestyle='--')
    plt.xlim([0.6,max_x])
    plt.ylim([-1,1]) 
    plt.xlabel('Ratio (observed/expected sites)')
    plt.tight_layout()
    plt.savefig('{}_motif_box_by_species_{}.pdf'.format(seq_type,method),dpi=500,bbox_inches='tight',transparent=True)
    plt.close()

def plot_density_curves(ratio_dict,method):
    plt.figure(figsize=(6,2.5))
    sns.set_style('white')
    colours = ['#f7b2b7','#f7717d','#de639a','#7f2982','#16001e','#92b4f4','#5e7ce2','#0a369d','#f7b05b','#af9b46']
    for i,rm_type in enumerate(sorted(ratio_dict.keys())):
        mean_per_species_motif = []
        for motif in ratio_dict[rm_type]:
            mean_per_species_motif.append(np.mean(ratio_dict[rm_type][motif]))
        print '{}: {} collapsed motifs'.format(rm_type,len(mean_per_species_motif))
        sns.distplot(mean_per_species_motif,label='{}, N = {}'.format(rm_type,len(mean_per_species_motif)),color=colours[i],hist=False,rug=True)
        #else:
        #    print '{} excluded, < 20 examples'.format(rm_type)
        #sns.distplot(ratio_dict[rm_type],label=rm_type,color=colours[i],hist=False)
    plt.plot([1,1],[0,50],linestyle='--',color='#D3D3D3')
    plt.xlim([0,3])
    plt.ylim([0,3])
    plt.xlabel('Ratio (observed/expected sites)')
    plt.ylabel('Density')
    plt.tight_layout()
    plt.legend()
    plt.savefig('mean_motif_distribution_by_type_{}.pdf'.format(method),dpi=500,bbox_inches='tight')
    plt.close()

    """#plot results for species
    plt.figure(figsize=(3,3))
    sns.set_style('white')
    plt.violinplot([observed_hits,markov_hits],showmeans=True, showextrema=True, showmedians=True)
    plt.xticks(range(1,3),[seq_type,'markov'])
    plt.ylabel('#/genome')
    plt.tight_layout()
    plt.savefig(seq_type+'_motifs_markov_'+bacteria[0] + bacteria.split()[1]+ '.pdf',dpi=500,bbox_inches='tight')
    plt.close()

    print 'motif enrichment:', sum(observed_hits),seq_type+' hits,',sum(shuffled_motif_hits),'shuffled hits,',wilcoxon(observed_hits,shuffled_motif_hits)
    plt.figure(figsize=(3,3))
    sns.set_style('white')
    plt.violinplot([observed_hits,shuffled_motif_hits],showmeans=True, showextrema=True, showmedians=True)
    plt.xticks(range(1,3),[seq_type,'shuffled'])
    plt.ylabel('#/genome')
    plt.tight_layout()
    plt.savefig(seq_type+'_motifs_shuffled_'+bacteria[0] + bacteria.split()[1]+ '.pdf',dpi=500,bbox_inches='tight')
    plt.close()"""

                #kr = 1
                #for i in range(len(mtase_target)):
                #    kr = kr*np.exp()
                #shuff = []
                #for n in range(10):
                #    shuff.append(len(re.findall(mtase_target,''.join(random.sample(sequence,len(sequence))))))
                #if len(re.findall(mtase_target[1:-1],sequence)) > 0 : 
                #    shuff.append(len(re.findall(mtase_target[:-1],sequence))*len(re.findall(mtase_target[1:],sequence))/len(re.findall(mtase_target[1:-1],sequence)))
                #shuffled_motif_hits.append(np.mean(shuff)) #len(re.findall(mtase_target,''.join(random.sample(sequence,len(sequence))))))
