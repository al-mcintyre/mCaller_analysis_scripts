abm237@node147 motif_underrepresentation $ grep -v "Type III" rusinov_RM_species_pairs.tsv | grep "Type II" | cut -f 4 | awk '{ml+=length($1)}END{print ml/NR}'
5.31269
abm237@node147 motif_underrepresentation $ grep -v "Type III" rusinov_RM_species_pairs.tsv | grep "Type II" | cut -f 4 | awk -v min=10 '{if (length($1) < min) min=length($1)}END{print min}'
2
abm237@node147 motif_underrepresentation $ grep -v "Type III" rusinov_RM_species_pairs.tsv | grep "Type II" | cut -f 4 | awk -v max=0 '{if (length($1) > max) max=length($1)}END{print max}'
12
abm237@node147 motif_underrepresentation $ grep -v "Type III" rusinov_RM_species_pairs.tsv | grep -v "Type II" | grep -v "Type IV" | grep "Type I" | cut -f 4 | awk -v max=0 '{if (length($1) > max) max=length($1)}END{print max}'
16
abm237@node147 motif_underrepresentation $ grep -v "Type III" rusinov_RM_species_pairs.tsv | grep -v "Type II" | grep -v "Type IV" | grep "Type I" | cut -f 4 | awk -v min=10 '{if (length($1) < min) min=length($1)}END{print min}'
10
abm237@node147 motif_underrepresentation $ grep -v "Type III" rusinov_RM_species_pairs.tsv | grep -v "Type II" | grep -v "Type IV" | grep "Type I" | cut -f 4 | awk '{ml+=length($1)}END{print ml/NR}'
13.4395

