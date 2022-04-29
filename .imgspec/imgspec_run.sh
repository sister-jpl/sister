#!/bin/bash

imgspec_dir=$(cd "$(dirname "$0")" ; pwd -P)
pge_dir=$(dirname ${imgspec_dir})

source activate sister

mkdir output temp
input_file=$(ls input/*.*)

python ${pge_dir}/scripts/l1_preprocess.py $input_file output/ temp/ $1

cd output
out_dir=$(ls ./)
mv ${out_dir} ${out_dir}_l1p
tar -czvf ${out_dir}_l1p.tar.gz ${out_dir}_l1p
rm -r $out_dir
