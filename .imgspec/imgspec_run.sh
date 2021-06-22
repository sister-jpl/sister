#!/bin/bash

# Get directories and paths for scripts
imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
sister_dir=$(dirname $imgspec_dir)

echo "imgspec_dir is $imgspec_dir"
echo "sister_dir is $isofit_dir"

# input/output dirs
input="input"
output="output"
mkdir -p ${output}
