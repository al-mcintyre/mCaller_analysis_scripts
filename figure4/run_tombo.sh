newdir=alb4_barcode03_copy
ref=../Efaecalis/efaecalis_pb.fasta
#source activate tomboenv
#tombo resquiggle $newdir $ref --processes 32

tombo test_significance --fast5-basedirs $newdir --dna \
    --alternate-bases 5mC --statistics-file-basename sample.m5C
tombo write_wiggles --wiggle-types fraction --statistics-filename sample.m5C.5mC.tombo.stats
mv tombo_results.fraction_modified_reads.minus.wig tombo_results.fraction_m5C_reads.minus.wig
mv tombo_results.fraction_modified_reads.plus.wig tombo_results.fraction_m5C_reads.plus.wig

tombo test_significance --fast5-basedirs $newdir \
    --alternate-bases 6mA --statistics-file-basename sample.m6A
tombo write_wiggles --wiggle-types fraction --statistics-filename sample.m6A.6mA.tombo.stats
mv tombo_results.fraction_modified_reads.minus.wig tombo_results.fraction_m6A_reads.minus.wig
mv tombo_results.fraction_modified_reads.plus.wig tombo_results.fraction_m6A_reads.plus.wig

