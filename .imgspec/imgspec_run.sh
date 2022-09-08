#!/bin/bash

imgspec_dir=$(cd "$(dirname "$0")" ; pwd -P)
pge_dir=$(dirname ${imgspec_dir})

source activate sister
mkdir output temp

input_file=input/*.*
input_file=$(ls input/*.* | tail -n 1)

echo $input_file
base=$(basename $input_file)

if [[ $base == PRS* ]]; then
    echo $2
    aws s3 cp $2 ./input
    lst_archive=$(ls input/*landsat.tar.gz)
    tar -xzvf $lst_archive -C input/
    landsat=$(ls input/*landsat)
    rdn_coeffs=${pge_dir}/data/prisma/*_radcoeff_surface.npz
    smile=${pge_dir}/data/prisma/*_wavelength_shift_surface_smooth.npz
    prisma_zip=input/*.zip
    python ${pge_dir}/scripts/l1_preprocess.py $prisma_zip output/ temp/ $1 $smile $rdn_coeffs $landsat
    datetime=$(echo $(sed 's/.\{6\}$/T&/' <<< $(echo $base | cut -c17-30)))
    file_base=SISTER_PRISMA_${datetime}_L1B
    rm output/*/*.log
    rm output/*/*.csv

else
    python ${pge_dir}/scripts/l1_preprocess.py $input_file output/ temp/ $1

    if [[ $base == DESIS* ]]; then
        datetime=$(echo $base | cut -c32-46)
        file_base=SISTER_DESIS_${datetime}_L1B
    fi
fi

# Rename PRISMA and DESIS files
if [[ $base == PRS* ]] || [[ $base == DESIS* ]]; then

    mv output/*/ output/${file_base}_RDN_000

    mv output/${file_base}_RDN_000/*rdn_prj.hdr output/${file_base}_RDN_000/${file_base}_RDN_000.hdr
    mv output/${file_base}_RDN_000/*rdn_prj output/${file_base}_RDN_000/${file_base}_RDN_000

    mv output/${file_base}_RDN_000/*obs_prj.hdr output/${file_base}_RDN_000/${file_base}_OBS_000.hdr
    mv output/${file_base}_RDN_000/*obs_prj output/${file_base}_RDN_000/${file_base}_OBS_000

    mv output/${file_base}_RDN_000/*loc_prj.hdr output/${file_base}_RDN_000/${file_base}_LOC_000.hdr
    mv output/${file_base}_RDN_000/*loc_prj output/${file_base}_RDN_000/${file_base}_LOC_000
fi

cd output

# Future....replace placeholder with CRID, bad practice runs twice first the foldername is changed then the files...
# should only be needed for AVIRIS data, CRID can be used when renaming DESIS and PRISMA imagery
# maybe not needed can pass CRID pge script and use for AVIRIS renaming
#find . -iname "*_000*" | rename 's/\_000/\_CRID/g';
#find . -iname "*_000*" | rename 's/\_000/\_CRID/g';

out_dir=$(ls ./)

tar -czvf ${out_dir}.tar.gz ${out_dir}
rm -r ${out_dir}
