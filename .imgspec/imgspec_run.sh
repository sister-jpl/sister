#!/bin/bash
# arg1: Output resolution in meters, should be a multiple of 30 (ex. 30,60,90.....)

imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
sister_dir=$(dirname ${imgspec_dir})

aws_cop_url='https://copernicus-dem-30m.s3.amazonaws.com/'
shift_surface='https://github.com/EnSpec/sister/raw/master/data/prisma/wavelength_shift/PRISMA_20200721104249_20200721104253_0001_wavelength_shift_surface'
output_dir='output'
temp_dir='tmp'
l1_zip=input/PRS*.zip

mkdir -p $output_dir
mkdir -p $temp_dir

# Run PRISMA PGE, export rdn, obs and loc ENVI files
python ${sister_dir}/scripts/prisma/prisma_pge.py $l1_zip $output_dir $temp_dir $aws_cop_url -proj -res $1 -shift $shift_surface

# gzip output files in preparation for downstream processing
cd $output_dir
l1_output_dir=PRS*
tar -cf ${l1_output_dir}.tar ${l1_output_dir}
rm -rf $l1_output_dir
gzip ${l1_output_dir}.tar