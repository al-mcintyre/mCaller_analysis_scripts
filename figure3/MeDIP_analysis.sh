table=species_motif_table.txt
uniq_species=$(cut -f 2 $table | sort | uniq)
#mkdir bwa_inds

for Species in Cneoformans; do #$uniq_species; do
    echo $Species
    Sp=$(grep $Species $table | head -1 | cut -f1)
    species=$(grep $Species $table | head -1 | cut -f3)
    genome_size=$(grep $Species $table | head -1 | cut -f6)
    ref=/athena/masonlab/scratch/projects/microbes/zymo_strains/assemblies/${species}_pb.fasta
    end=pb
    if [ ! -f $ref ]; then 
        ref=/athena/masonlab/scratch/projects/microbes/zymo_strains/assemblies/${species}_ontlumina.fasta
        end=ontlumina
    fi 
    if [ ! -f $ref ]; then 
        ref=/athena/masonlab/scratch/projects/microbes/zymo_strains/assemblies/${species}_ont.fasta
        end=ont
    fi
    
    #bwa index $ref -p bwa_inds/${species}_${end}
    dir=$Species
    #mkdir dir; mv ${Sp}*.fasta.gz $dir

    #bwa mem -t 8 bwa_inds/${species}_${end} $dir/${Sp}_input.fastq.gz | samtools view -bS > $dir/${Sp}_input.bam
    #samtools sort -T ${Sp}_input_tmp.bam -o $dir/${Sp}_input.sorted.bam $dir/${Sp}_input.bam
    #samtools index $dir/${Sp}_input.sorted.bam
    #rm $dir/${Sp}_input.bam
    #bwa mem -t 8 bwa_inds/${species}_${end} $dir/${Sp}_IP.fastq.gz | samtools view -bS > $dir/${Sp}_IP.bam
    #samtools sort -T ${Sp}_IP_tmp.bam -o $dir/${Sp}_IP.sorted.bam $dir/${Sp}_IP.bam
    #samtools index $dir/${Sp}_IP.sorted.bam
    #rm $dir/${Sp}_IP.bam

    #macs2 callpeak -t $dir/${Sp}_IP.sorted.bam -c $dir/${Sp}_input.sorted.bam --nomodel --extsize 100 -g $genome_size -n $dir/${Sp}.macs2 -f BAM --verbose 3

    motifs=$(grep $Species $table | cut -f 4)
    posfi=../../../20171020_zymo_strains_tmp/${Species}/${species}_positions_to_check.txt
    ontfi=${Species}_ont.bed
    if [ -f ../../../20171020_zymo_strains_tmp/${Species}/${species}.methylation.summary.bed ]; then
    bedtools sort -i ../../../20171020_zymo_strains_tmp/${Species}/${species}.methylation.summary.bed > $ontfi

    for motif in $motifs; do
        if [ "$motif" != "none" ]; then
            for status in canon missed; do 
                cut -f 1-6,8-9 $posfi | grep -v noncanon | grep $status | grep $motif | awk -v OFS="\t" '{print $1,$2,$3,".",".",$4}' | bedtools sort -i - > ${Species}/motif_positions.bed
                motifsoverlappingpeaks=$(bedtools intersect -a ${Species}/${Sp}.macs2_peaks.narrowPeak -b ${Species}/motif_positions.bed -wb | wc -l | cut -d' ' -f 1)
                motifsoverlappingont=$(bedtools intersect -a $ontfi -b ${Species}/motif_positions.bed -wb | wc -l | cut -d' ' -f 1)
                if [ "$status" = "canon" ]; then
                    fullinter=$(bedtools multiinter -i $ontfi ${Species}/motif_positions.bed ${Species}/${Sp}.macs2_peaks.narrowPeak | awk '($4 == 3){print $0}' | wc -l | cut -d' ' -f 1)
                else 
                    fullinter=$(bedtools multiinter -i ${Species}/motif_positions.bed $ontfi ${Species}/${Sp}.macs2_peaks.narrowPeak | awk '($4 == 1 && $5 == 1){print $0}' | wc -l | cut -d' ' -f 1)
                fi
                totmotifs=$(wc -l ${Species}/motif_positions.bed | cut -d' ' -f 1)
                if [ "$totmotifs" -gt 0 ]; then
                    echo " - " $motif $status ': ' $motifsoverlappingpeaks 'out of' $totmotifs 'motifs found within peaks,' $motifsoverlappingont 'ONT,' $fullinter 'from full intersect'
                fi
            done
        fi
        cut -f 1-6,8-9 $posfi | grep -v control | grep $motif | awk -v OFS="\t" '{print $1,$2,$3,".",".",$4}' | bedtools sort -i - > ${Species}/motif_positions.bed
        motifsoverlappingpeaks=$(bedtools intersect -a ${Species}/${Sp}.macs2_peaks.narrowPeak -b ${Species}/motif_positions.bed -wb | wc -l | cut -d' ' -f 1)
        motifsoverlappingont=$(bedtools intersect -a $ontfi -b ${Species}/motif_positions.bed -wb | wc -l | cut -d' ' -f 1)
        fullinter=$(bedtools multiinter -i $ontfi ${Species}/motif_positions.bed ${Species}/${Sp}.macs2_peaks.narrowPeak | awk '($4 == 3){print $0}' | wc -l | cut -d' ' -f 1)
        totmotifs=$(wc -l ${Species}/motif_positions.bed | cut -d' ' -f 1)
        echo " - " $motif 'all: ' $motifsoverlappingpeaks 'out of' $totmotifs 'motifs found within peaks,' $motifsoverlappingont 'ONT,' $fullinter 'from full intersect'
    done
    grep -v noncanon ../../../20171020_zymo_strains_tmp/${Species}/${species}_positions_to_check.txt | grep -v control | grep -v species | grep -v none | awk -v OFS="\t" '{print $1,$2,$3,".",".",$4}' > ${Species}/motif_positions.bed
    peaksoverlappingmotifs=$(bedtools intersect -a ${Species}/${Sp}.macs2_peaks.narrowPeak -b ${Species}/motif_positions.bed -wa | cut -f 4 | uniq | wc -l | cut -d' ' -f 1)
    totpeaks=$(wc -l ${Species}/${Sp}.macs2_peaks.narrowPeak | cut -d' ' -f 1)
    echo " - " $peaksoverlappingmotifs 'out of' $totpeaks 'peaks overlap known motifs (missed and found by PacBio)'
    fi

    if [ -f ../../../20171207_2_zymo/${Species}/${species}.methylation.summary.bed ] || [ "$Species" == "Scerevisiae" ] ; then
        if [ -f ../../../20171020_zymo_strains_tmp/${Species}/${species}.methylation.summary.bed ]; then
            bedtools sort -i ../../../20171020_zymo_strains_tmp/${Species}/${species}.methylation.summary.bed > $ontfi
        else
            bedtools sort -i ../../../20171207_2_zymo/${Species}/${species}.methylation.summary.bed > $ontfi
        fi
        bedtools intersect -a $ontfi -b ${Species}/${Sp}.macs2_peaks.narrowPeak -wa | awk -v OFS="\t" '{print $1,$2-20,$3+20}' | bedtools getfasta -fi $ref -bed - > ${Species}_medip_ont_overlap.fasta
        bedtools intersect -a $ontfi -b ${Species}/${Sp}.macs2_peaks.narrowPeak -wb | cut -f 1,2,3 | sort | uniq | wc -l 

        #macs2 callpeak -t $dir/${Sp}_input.sorted.bam -c $dir/${Sp}_IP.sorted.bam --nomodel --extsize 100 -g $genome_size -n $dir/TN_${Sp}.macs2 -f BAM --verbose 3
        bedtools intersect -a $dir/TN_${Sp}.macs2_peaks.narrowPeak -b $ontfi -wa -v | awk -v OFS="\t" '{print $1,int(($2+$3)/2)-20,int(($2+$3)/2)+21}' | bedtools getfasta -fi $ref -bed - > ${Species}_medip_ont_control.fasta 

        head ${Species}_medip_ont_overlap.fasta
        head ${Species}_medip_ont_control.fasta

        mkdir ${Species}_homer_motifs
        findMotifs.pl ${Species}_medip_ont_overlap.fasta fasta ${Species}_homer_motifs #-fasta ${Species}_medip_ont_control.fasta

    fi

done
