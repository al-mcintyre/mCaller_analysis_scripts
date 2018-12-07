#!/bin/bash

#this script generates mCaller output files for all quality levels and with max 2 skips -- these can be further processed during plotting

#usage: mCaller_nanopolish.py [-h] (-p POSITIONS | -m) -r REFERENCE -f TSV
#                             [-t THREADS] [-l LABEL] [-b BASE] [--train] [-v]
#REF=/athena/masonlab/scratch/projects/nasa/biomolecule_sequencer/data/pacbio/pb_ecoli_polished_assembly.fasta
REF=/pbtech_mounts/masonlab_scratch_athena/projects/nasa/biomolecule_sequencer/data/pacbio/pb_ecoli_polished_assembly.fasta
TSV='/athena/masonlab/scratch/users/abm237/mCaller_old_analyses/mCaller_bioRxiv_scripts_figures/chiu_R9_fwd.gm.ecoli.eventalign.tsv'
#TSV='/athena/masonlab/scratch/users/abm237/mCaller_old_analyses/R9_nasa_data_chiu/chiu_R9.bwa.ecoli.eventalign.tsv' #need to change read name format 
#TESTSV='/athena/masonlab/scratch/users/abm237/mCaller_old_analyses/mCaller_bioRxiv_scripts_figures/mason_R9_fwd.gm.ecoli.eventalign.tsv'
TESTSV=m6A_test/mason_R9_fwd.gm.ecoli.eventalign.tsv #moved from previous location
POS=ecoli_training_list.txt
PO='ecoli_training_'
FASTQ='/athena/masonlab/scratch/users/abm237/mCaller_old_analyses/R9_nasa_data_chiu/chiu_R9_fwd.fastq'
TASTQ='/athena/masonlab/scratch/users/abm237/mCaller_old_analyses/R9_nasa_data_mason/20160922_nasa_R9_fwd.fastq'
MODFI='m6A_model_'
TT=chiu_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.train #.reformat #'/athena/masonlab/scratch/users/abm237/mCaller_old_analyses/mCaller_bioRxiv_scripts_figures/chiu_R9_fwd.gm.ecoli.eventalign.diffs'
TOS='ecoli_testing_'
#POS=second_half.txt
#TSV='mason_R9_fwd.gm.ecoli.eventalign.tsv'

LEN=6
C=NN
## for comparison to old results: awk -v OFS="\t" '{print "ecoli",$2,$5,$6}' ../mCaller_bioRxiv_scripts_figures/chiu_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.[mA]* | sort | uniq > old_training_set.txt
#python make_training_sets.py 
#cut -f 1,2,4,6 ecoli_training.txt > $POS
#cut -f 1,2,4,6 ecoli_training.txt | grep m6A > ${PO}m6A.txt
#cut -f 1,2,4,6 ecoli_training.txt | grep -v m6A > ${PO}A.txt
#cut -f 1,2,4,6 ecoli_testing.txt | grep m6A > ${TOS}m6A.txt
#cut -f 1,2,4,6 ecoli_testing.txt | grep -v m6A > ${TOS}A.txt

#have to set twobase = False in mCaller extraction file 
:<<'END'
for LEN in 4 8; do 
    #time ~/bin/mCaller/mCaller.py -p $POS -r $REF -e $TSV -t 8 -b A -n $LEN --train -f $FASTQ -d ${MODFI}_${C}_${LEN}_wskips.pkl -c $C -s 0 #increase skips to 2 for LEN=6

    time ~/bin/mCaller/mCaller.py -p ${TOS}A.txt -r $REF -e $TESTSV -t 8 -b A -n $LEN -f $TASTQ -d ${MODFI}_${C}_${LEN}_wskips.pkl #-s 2
    mv m6A_test/mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN} mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.A
    time ~/bin/mCaller/mCaller.py -p m6A_test/${TOS}m6A.txt -r $REF -e $TESTSV -t 8 -b A -n $LEN -f $TASTQ -d ${MODFI}_${C}_${LEN}_wskips.pkl #-s 2
    mv m6A_test/mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN} mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.m6A
    awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,"m6A"}' mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.m6A > mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.m6A2; mv mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.m6A2 mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.m6A
    awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,"A"}' mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.A > mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.A2; mv mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.A2 mason_R9_fwd.gm.ecoli.eventalign.diffs.${LEN}.A

done 
END

#time NN vs. RF (model_7_NN_chiu.pkl vs model_7_RF_chiu.pkl)
for MOD in RF NN; do
    echo $MOD
    time ~/bin/mCaller/mCaller.py -p ${TOS}A.txt -r $REF -e $TESTSV -t 8 -b A -n $LEN -f $TASTQ -d model_7_${MOD}_chiu.pkl
done

