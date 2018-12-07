threads=6

ref=bsubtilis_pb.fasta 
filename=bsubtilis
bar=05
echo $ref", "$filename
if [ ! -f $filename.fastq ]; then
    for dir in pass; do 
        d=../albacore4/workspace/$dir/barcode$bar/
        cat $d/*.fastq* | bioawk -c fastx '{if (meanqual($qual) > 10) print "@"$name"\t"$comment"\n"$seq"\n+\n"$qual}'  > ${filename}.fastq
        python ~/SCRIPTS/remove_duplicates.py ${filename}.fastq
        mv ${filename}_clean.fastq ${filename}.fastq
        nanopolish index -d $d ${filename}.fastq
    done
fi

if [ ! -f $filename.sorted.bam ]; then
    graphmap align -t $threads -r $ref -d ${filename}.fastq -o $filename.sam 
    echo 'aligned'
    samtools view -bS $filename.sam | samtools sort -T /tmp/$filename.sorted -o $filename.sorted.bam 
    rm $filename.sam
    samtools index $filename.sorted.bam 
fi
if [ ! -f $filename.eventalign.tsv ]; then 
    nanopolish eventalign -t $threads --scale-events -n -r ${filename}.fastq -b $filename.sorted.bam -g $ref > $filename.eventalign.tsv 
    echo 'nanopolish eventaligned'
fi
model='/home/abm237/bin/mCaller/r95_twobase_model_NN_6_m6A.pkl'
if [ ! -f $filename.eventalign.diffs.6 ]; then
    /usr/bin/time /home/abm237/bin/mCaller/mCaller.py -t $threads -m A -r $ref -e $filename.eventalign.tsv -d $model -f $filename.fastq
    echo 'mCaller done'
fi

/home/abm237/bin/mCaller/make_bed.py -d 15 -t 0.5 -f $filename.eventalign.diffs.6
/home/abm237/bin/mCaller/make_bed.py -d 15 -t 0.5 -f $filename.eventalign.diffs.6 --control
echo 'beds made'

sort -rk1,1n -rk2,2n $filename.methylation.summary.bed | cut -f 1,2,3,4,6 | awk '{id+=1; print ">"$1":"$2-5":"$3+4":"$5":"id"\n"$4}' | sed 's/M/A/g' > $filename.methylated.fasta
sort -rk1,1n -rk2,2n $filename.methylation.control.summary.bed | cut -f 1,2,3,4,6 | shuf -n 100000 | awk '{id+=1; print ">"$1":"$2-5":"$3+4":"$5":"id"\n"$4}' | sed 's/M/A/g' > $filename.methylated.control.fasta

for exp in methylation methylation.control; do
    awk -v OFS="\t" '{print $1,$2-20,$3+20,".",".",$6}' $filename.$exp.summary.bed  | shuf -n 100000 > test.bed;  bedtools getfasta -s -fi $ref -bed test.bed > $filename.$exp.expanded.fasta
done
ame --o ame_all --control ${filename}.methylation.control.expanded.fasta ${filename}.methylation.expanded.fasta ../../REFS/rebase/rebase_all_M.meme
