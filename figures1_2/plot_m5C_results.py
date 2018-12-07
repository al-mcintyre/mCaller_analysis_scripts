from classifier import *
import sys
import random
import cPickle
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from decimal import Decimal
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from scipy.stats.stats import pearsonr # use spearman instead of pearson
from scipy.stats import spearmanr 

def extract_data(exp,num,i,qual_thresh,skip_thresh,mod_thresh,diffs_by_context,diff_fis=None,update=False):
   signals, contexts, labels, positions = [],[],[],[]
   if not diff_fis:
      diffs_fis = [exp+'_R9_fwd.gm.ecoli.eventalign.diffs.'+num+'m6A',exp+'_R9_fwd.gm.ecoli.eventalign.diffs.'+num+'A']
   for diff_fi in diff_fis:
      for line in open(diff_fi,'r'):
         chrom, readname, pos, context, diffs, strand, label = line.split('\t')
         read_qual = float(diffs.split(',')[-1])
         try:
          if read_qual >= qual_thresh:
            nskips = len([float(x) for x in diffs.split(',') if float(x) == 0])
            nmods = len([context.split('M')]) -1
            features = [float(x) for x in diffs.split(',')]
            if len(features) == i and nskips <= skip_thresh and nmods <= mod_thresh:
               if context not in diffs_by_context[label.strip()]:
                  diffs_by_context[label.strip()][context] = [features]
               else:
                  diffs_by_context[label.strip()][context].append(features)
               signals.append(features)
               contexts.append(context)
               labels.append(label.strip())
               positions.append(pos)
         except IndexError:
           continue
   for label in diffs_by_context:
      print len(diffs_by_context[label]),'contexts for',label
   print len(set(diffs_by_context[unmeth].keys())&set(diffs_by_context[meth].keys())),'overlapping contexts'
   print len(signals), 'observations' #, eg.',training_signals[:5]
   if update:
      return diffs_by_context, signals, contexts, labels, positions
   else:
      return signals, contexts, labels, positions

def plot_diffs_by_context(diffs_by_context,base,methbase,cols):
   count = 0
   sns.set_style('white')
   for context in diffs_by_context[methbase]:
      if context in diffs_by_context[base]:
         fig = plt.figure(figsize=(6,4))
         count += 1
         for label in diffs_by_context:
            col = cols[label]
            for i,entry in enumerate(diffs_by_context[label][context]):
               if i == 0:
                  plt.plot(range(1,len(entry[:-1])+1),entry[:-1],color=col,alpha=1,label=label)
               else:
                  plt.plot(range(1,len(entry[:-1])+1),entry[:-1],color=col,alpha=1)
         plt.title(context)
         plt.ylim([-8,8])
         sns.despine(trim=False)
         plt.ylabel('Measured - Expected Current (pA)')
         plt.xlabel('Cytosine Position in 6-mer')
         plt.legend()
         plt.savefig('R9_measured_vs_expected_across_positions_'+context+'.pdf',dpi=500,bbox_inches='tight',transparent=True)
      if count > 10:
         break

   fig = plt.figure(figsize=(6,4))
   for lab in [meth,unmeth]:
      col = cols[lab]
      all_contexts = diffs_by_context[lab].keys()
      contexts = random.sample(all_contexts,500)
      for i,context in enumerate(contexts):
        for entry in random.sample(diffs_by_context[lab][context],1):
          if i == 0:
             plt.plot(range(1,len(entry[:-1])+1),entry[:-1],color=col,alpha=0.2,label=lab)
          else:
             plt.plot(range(1,len(entry[:-1])+1),entry[:-1],color=col,alpha=0.2)
   plt.ylabel('Measured - Expected Current (pA)')
   plt.xlabel('Cytosine Position in 6-mer')
   plt.legend()
   sns.despine(trim=False)
   plt.ylim([-15,15])
   plt.savefig('R9_measured_vs_expected_across_positions.pdf',transparent=True,dpi=500,bbox_inches='tight')

   fig = plt.figure(figsize=(6,4))
   changes_by_position = {'position':[],'base':[],'diff':[]}
   for lab in diffs_by_context:
      for context in diffs_by_context[lab]:
         for entry in diffs_by_context[lab][context]:
            for pos,diff in enumerate(entry[:-1]):
               changes_by_position['position'].append(pos+1)
               changes_by_position['base'].append(lab)
               changes_by_position['diff'].append(diff)
   dPos = pd.DataFrame(changes_by_position)
   sns.boxplot(x="position", y="diff", hue="base", data=dPos, palette=[cols[base],cols[methbase]])
   sns.despine(trim=False)
   plt.xlabel('Cytosine Position in 6-mer')
   plt.ylabel('Measured - Expected Current (pA)')
   plt.ylim([-40,40])
   plt.legend(title='',loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True)
   plt.savefig('R9_change_by_position_box.png',transparent=True,dpi=500, bbox_inches='tight')

unmeth = 'C'
meth = 'm5C'
roc = {}
prob_scores = {}
nums = ['6.']
num_feats = [7]
classifiers = ['NN']
base_colours = {unmeth:'#55B196', meth:'#B4656F'}
diffs_by_context = {unmeth:{},meth:{}}
training_signals,training_labels,training_positions,training_contexts = defaultdict(list),defaultdict(list),defaultdict(list),defaultdict(list)
testing_signals,testing_labels,testing_positions,testing_contexts = defaultdict(list),defaultdict(list),defaultdict(list),defaultdict(list)
saveas='R9_m5C'

for exp,modification in zip(['mssl','pcr'],[meth,unmeth]):
   diff_fis= ['timp_ecoli_'+exp+'.eventalign.diffs.'+nums[0]+modification+'2']
   
   #for qual, qual_threshold in zip(['hq.',''],[9,0]):
   for num,i in zip(nums,num_feats):
      if num == '6.':
         quals,qthreshs = [''],[0] # ['_mq',''],[7,0]
         skips,sthreshs = [''],[0] #['_sk',''],[2,0]
         mods,modthreshs = [''],[100] #['_1mod',''],[1,i*2-1] #interestingly, no ROC plot difference in predictions for 11mers with multiple m5Cs and single but I don't think the threshold worked
      else:
         quals,qthreshs = [''],[0]
         skips,sthreshs = [''],[0]
      for qual, qual_threshold in zip(quals,qthreshs):
        for mod, mod_threshold in zip(mods,modthreshs):
         for skip, skip_threshold in zip(skips,sthreshs):
            clf = 'NN'+qual+mod
            print exp, num, clf
            if quals == '':
               diffs_by_context,signals,contexts,labels,positions = extract_data(exp,num,i,qual_threshold,skip_threshold,mod_threshold,diffs_by_context,diff_fis,update=True)
            else:
               signals,contexts,labels,positions = extract_data(exp,num,i,qual_threshold,skip_threshold,mod_threshold,diffs_by_context,diff_fis)
      
            data_split = len(signals)/10

            training_signals[clf] = training_signals[clf] + signals[:data_split]
            testing_signals[clf] = testing_signals[clf] + signals[data_split:]
            training_labels[clf] = training_labels[clf] + labels[:data_split]
            testing_labels[clf] = testing_labels[clf] + labels[data_split:]
            training_positions[clf] = training_positions[clf] + positions[:data_split]
            testing_positions[clf] = testing_positions[clf] + positions[data_split:]
            training_contexts[clf] = training_contexts[clf] + contexts[:data_split]
            testing_contexts[clf] = testing_contexts[clf] + contexts[data_split:]

plot_diffs_by_context(diffs_by_context,unmeth,meth,base_colours)
sys.exit(0)
for clf in training_signals:
   for num in nums:
      i = len(training_signals[clf][0])
      print 'training model...',clf
      model = model_signal(training_signals[clf],training_labels[clf],True,modelfile='model_'+str(i)+'_'+clf+'.pkl',classifier=clf)  
      sys.stdout.flush()
      if clf not in prob_scores:
         prob_scores[clf] = {}
      print 'classifying test sites...'
      prob_scores[clf][num],results = model_signal(testing_signals[clf],testing_labels[clf],False,modelfile='model_'+str(i)+'_'+clf+'.pkl')

      accuracy = []
      for base in prob_scores[clf][num]:
         if base == unmeth:
            ind = 0 
         else:
            ind = 1
         accuracy.append( len([x for x in prob_scores[clf][num][base] if x[ind] > 0.5])*100./len(prob_scores[clf][num][base]) )
         print base, accuracy[-1]
      print np.mean(accuracy) 

      sns.set_style('darkgrid')
      sns.set_palette(['#55B196','#B4656F'])
      fig = plt.figure(figsize=(3,4))
      prob_dict = {'probability':[x[1] for x in prob_scores[clf][num][unmeth]]+[x[1] for x in prob_scores[clf][num][meth]],'base':[unmeth]*len(prob_scores[clf][num][unmeth])+[meth]*len(prob_scores[clf][num][meth])}
      prob_db = pd.DataFrame(prob_dict)
      print '\nplotting ' + saveas+'.'+clf+'.'+num+'png'
      sns.boxplot(x="base", y="probability", data=prob_db)
      sns.despine()
      plt.savefig(saveas+'.'+clf+'.'+num+'png',transparent=True,dpi=500,bbox_inches='tight')
      plt.close()

      if clf not in roc:
         roc[clf] = {}
      roc[clf][num] = {'tp':[],'fp':[]}

              #print prob_scores['A'][:10]
      for thresh in np.arange(0,1.1,0.1)[::-1]:
         roc[clf][num]['tp'].append(len([x for x in prob_scores[clf][num][meth] if x[1]>=thresh])*100./len(prob_scores[clf][num][meth]))
         roc[clf][num]['fp'].append(len([x for x in prob_scores[clf][num][unmeth] if x[1]>=thresh])*100./len(prob_scores[clf][num][unmeth]))
              #print roc[clf][num]['tp'], roc[clf][num]['fp'], np.max(prob_scores['m6A']), np.max(prob_scores['A']) 
         
      pos_dict = {unmeth:{},meth:{}}
      """if clf == 'NN' and num == '6.':
         roc[clf+'_pos'] = {}
         roc[clf+'_pos'][num] = {'tp':[],'fp':[]}
         for pos,lab,res in zip(positions,labels,results):
            if pos not in pos_dict[lab]:
               pos_dict[lab][pos] = []
            if res[0] > 0.5:
               read_score = 1
            else:
               read_score = 0
            pos_dict[lab][pos].append(read_score)
         for thresh in np.arange(0,1,0.1)[::-1]:
            roc[clf+'_pos'][num]['tp'].append(len([pos for pos in pos_dict[meth] if np.mean(pos_dict[meth][pos])<=thresh])*100./len(pos_dict[meth]))
            roc[clf+'_pos'][num]['fp'].append(len([pos for pos in pos_dict[unmeth] if np.mean(pos_dict[unmeth][pos])<=thresh])*100./len(pos_dict[unmeth]))
      """

sns.set_style('white')
#nums = ['6.','10.','15.']
num_lab ={'4.':'7','8.':'15','':'11','6.':'11','10.':'19','15.':'29'}
clf_colours = {'NN_mq':['#D15C9C','#D15C9C','#D15C9C'],'NN':['#D15C9C','#B70064','#750040'],'NN_1mod':['#750040','#750040','#750040'],'NN_mq_1mod':['#6C698D','#6C698D','#6C698D'],'SVM':['#E1ADE1','#C382C3','#754E75'],'NN_hq':['#E1ADE1','#C382C3','#754E75'],'RF':['#96B596','#5B8C5A','#3A5A3A'],'NBC':['#B9B9BB','#767677','#000009'],'LR':['#87BCDE','#6F9AB6','#3E5665'],'NN_sk':['#FF5C72','#FF0022','#A30016'],'NN_pos':['#A19FB6','#6C698D','#45435A']}

num = '6.'
num_ind = 1
plt.figure(figsize=(4,4))
for clf in training_signals: #['NN','NN_hq','NN_mq']:
   plt.plot(roc[clf][num]['fp'],roc[clf][num]['tp'],color=clf_colours[clf][num_ind],label=clf)
plt.legend(loc=4)
plt.tight_layout()
plt.ylabel('True Positives')
plt.xlabel('False Positives')
plt.ylim([0,100])
plt.xlim([0,100])
plt.savefig(saveas+'.ROC.pdf',dpi=500,transparent=True,bbox_inches='tight')
plt.close()
