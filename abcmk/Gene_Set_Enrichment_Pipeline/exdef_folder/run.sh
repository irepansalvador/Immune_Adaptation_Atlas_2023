#!/bin/bash

case=${1}

cp params.txt params/params_${case}
cp make_input.txt params/make_${case}

sed -i "s/case/${case}/" params/params_${case}
sed -i "s/case/${case}/" params/make_${case}

bash distance.sh -c datasets/${case}.txt -a ensembl_gene_coords_v109.bed
perl make_valid_file.pl params/make_${case} params/valid_${case}

#julia bootstrap_cmd.jl -v params/valid_${case} -d params/distance_${case} -m 1 -r 3 -c datasets/${case}  -f jesus_g3_table_ready.txt -o outputs/${case} -t 0.05 -i 10 -n 10 -p

perl exdef_pipeline.pl params/params_${case}

cat VIPs/file_1  >  outputs/${case}_case.txt
cat nonVIPs/file_1 > outputs/${case}_control.txt

rm params/params_${case} params/make_${case} params/valid_${case}  params/distance_${case}.txt
