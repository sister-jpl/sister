#!/bin/bash

imgspec_dir=$(cd "$(dirname "$0")" ; pwd -P)
pge_dir=$(dirname ${imgspec_dir})

source activate sister
mkdir output temp

shopt -s extglob
input_file=./input/!(*landsat*)

base=$(basename $input_file)

if [[ $input_file == PRS* ]]; then
    landsat==$(ls input/*landsat)
    rdn_coeffs=${pge_dir}/data/prisma/*_radcoeff_surface.npz
    smile=${pge_dir}/data/prisma/*_wavelength_shift_surface_smooth.npz
    python ${pge_dir}/scripts/l1_preprocess.py $input_file output/ temp/ $1 $smile $rdn_coeffs
else
    python ${pge_dir}/scripts/l1_preprocess.py $input_file output/ temp/ $1
fi

cd output
out_dir=$(ls ./)
mv ${out_dir} ${out_dir}_l1p
tar -czvf ${out_dir}_l1p.tar.gz ${out_dir}_l1p
rm -r ${out_dir}_l1p
