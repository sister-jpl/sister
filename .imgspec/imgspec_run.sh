#!/bin/bash
# arg1: PRISMA L1 zip file
# arg2: Output directory
# arg3: Temporary directory

imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
sister_dir=$(dirname ${imgspec_dir})

aws_cop_url='https://copernicus-dem-30m.s3.amazonaws.com/'
shift_surface='https://github.com/EnSpec/sister/raw/master/data/prisma/wavelength_shift/PRISMA_20200721104249_20200721104253_0001_wavelength_shift_surface'

python ${sister_dir}/scripts/prisma/prisma_pge.py $1 $2 $3 ${aws_cop_url} -proj -res 30 -shift ${shift_surface}
