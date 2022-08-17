#!/bin/bash

imgspec_dir=$(cd "$(dirname "$0")" ; pwd -P)
pge_dir=$(dirname ${imgspec_dir})

source activate sister
mkdir output temp

input_file=input/*.*
echo $input_file
base=$(basename $input_file)

if [[ $base == PRS* ]]; then
    aws s3 cp $2 ./input
    lst_archive=$(ls input/*landsat.tar.gz)
    tar -xzvf $lst_archive -C input/
    landsat=$(ls input/*landsat)
    rdn_coeffs=${pge_dir}/data/prisma/*_radcoeff_surface.npz
    smile=${pge_dir}/data/prisma/*_wavelength_shift_surface_smooth.npz
    prisma_zip=input/*.zip
    python ${pge_dir}/scripts/l1_preprocess.py $prisma_zip output/ temp/ $1 $smile $rdn_coeffs $landsat
else
    python ${pge_dir}/scripts/l1_preprocess.py $input_file output/ temp/ $1
fi

cd output
out_dir=$(ls ./)
mv ${out_dir} ${out_dir}_l1p
tar -czvf ${out_dir}_l1p.tar.gz ${out_dir}_l1p
rm -r ${out_dir}_l1p
