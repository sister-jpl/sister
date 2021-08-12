
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import os
import hytools_lite as htl
from hytools_lite.io.envi import WriteENVI
import numpy as np
from scipy.interpolate import interp1d
from skimage.util import view_as_blocks

def main():
    ''' Perform a two-step spectral resampling to 10nm. First wavelengths are aggregated
    to approximatley 10nm, then aggregated spectra a interpolated to exactly 10nm using a
    cubic piecewise interpolator.
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str)
    parser.add_argument('outdir', type=str)
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    out_image = args.outdir + '/' + os.path.basename(args.infile) + "_10nm"
    hy_obj = htl.HyTools()
    hy_obj.read_file(args.infile,'envi')

    if hy_obj.wavelengths.max()< 1100:
        new_waves = np.arange(400,991,10)
    else:
        new_waves = np.arange(400,2451,10)

    bins = int(np.round(10/np.diff(hy_obj.wavelengths).mean()))
    agg_waves  = np.nanmean(view_as_blocks(hy_obj.wavelengths[:(hy_obj.bands//bins) * bins],
                                           (bins,)),axis=1)
    if args.verbose:
        print("Aggregating every: %s" % bins)

    out_header = hy_obj.get_header()
    out_header['bands'] = len(new_waves)
    out_header['wavelength'] = new_waves.tolist()
    out_header['fwhm'] = [10 for x in new_waves]
    out_header['default bands'] = []

    writer = WriteENVI(out_image,out_header)
    iterator =hy_obj.iterate(by = 'line')

    while not iterator.complete:
        if (iterator.current_line%100 == 0) and args.verbose:
            print(iterator.current_line)
        line = iterator.read_next()[:,:(hy_obj.bands//bins) * bins]
        line  = np.nanmean(view_as_blocks(line,(1,bins,)),axis=(2,3))
        interpolator = interp1d(agg_waves,line,fill_value = 'extrapolate', kind = 'cubic')
        line = interpolator(new_waves)
        writer.write_line(line,iterator.current_line)


if __name__ == "__main__":
    main()
