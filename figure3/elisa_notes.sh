#get total adenines per assembly
for fi in ../assemblies/*; do echo $fi; awk -v RS="A" 'END{print NR-1}'  $fi; awk -v RS="T" 'END{print NR-1}'  $fi; done 
