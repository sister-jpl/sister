#!/bin/bash

imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
sister_dir=$(dirname ${imgspec_dir})

aws_cop_url='https://copernicus-dem-30m.s3.amazonaws.com/'
shift_surface='https://github.com/EnSpec/sister/raw/master/data/prisma/wavelength_shift/PRISMA_20200721104249_20200721104253_0001_wavelength_shift_surface'
output_dir='/data2/prisma/rdn/'
temp_dir='/data2/temp/'
l1_zip=${sister_dir}/input/PRS*.zip

# Run PRISMA PGE, export rdn, obs and loc ENVI files
python ${sister_dir}/scripts/prisma/prisma_pge.py ${l1_zip} ${output_dir} ${temp_dir} ${aws_cop_url} -proj -res 30 -shift ${shift_surface}
