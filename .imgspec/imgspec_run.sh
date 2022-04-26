#!/bin/bash

imgspec_dir=$(cd "$(dirname "$0")" ; pwd -P)
pge_dir=$(dirname ${imgspec_dir})

conda create -y --name sister python=3.8
source activate sister
conda install -yc conda-forge gdal
conda install -y numba

cd $pge_dir
python setup.py install
yes | pip uninstall ray
yes | pip install 'ray[default]'

cd $imgspec_dir

mkdir output
mkdir temp

input_file=$(ls input/*.*)

python ${pge_dir}/scripts/l1_preprocess.py $input_file output/ temp/

cd output
out_dir=$(ls ./)
tar -czvf $out_dir.tar.gz $out_dir
rm -r $out_dir
