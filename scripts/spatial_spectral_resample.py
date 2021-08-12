import hytools_lite as htl
from skimage.util import view_as_blocks
import matplotlib.pyplot as plt
import numpy as np
from hytools_lite.io.envi import WriteENVI,envi_header_dict
from scipy.interpolate import interp1d


in_image = '/data2/prism/prm20160916t225841_rdn_v1w2/prm20160916t225841_rdn_v1w2_img'
hy_obj = htl.HyTools()
hy_obj.read_file(in_image,'envi')
block_size = 4


out_header = hy_obj.get_header()
out_header['lines']= hy_obj.lines//block_size
out_header['samples']=hy_obj.columns//block_size
out_header['data ignore value'] = -9999
# out_header['map info'] = map_info

spatial = '/data2/temp/prism_spatial_resample'
writer = WriteENVI(spatial,out_header)
iterator =hy_obj.iterate(by = 'band')

while not iterator.complete:
    band = np.copy(iterator.read_next())[:4* (img.lines//block_size),:4* (img.columns//block_size)]
    band[band == img.no_data] = np.nan
    band  = np.nanmean(view_as_blocks(band, (block_size,block_size)),axis=(2,3))
    band[np.isnan(band)] = -9999
    writer.write_band(band,iterator.current_band)





in_image = '/data2/desis/rfl/DESIS_DT0389141848_017-20191125T035518-V0213/DESIS_DT0389141848_017-20191125T035518-V0213_rfl_prj'
hy_obj = htl.HyTools()
hy_obj.read_file(in_image,'envi')


new_waves = np.arange(400,991,10)
bin_dict= {"DESIS": 3}
sensor = "DESIS"
agg_waves  = np.nanmean(view_as_blocks(hy_obj.wavelengths,(bin_dict[sensor],)),axis=1)

spectral = '/data2/temp/desis_spectral_test'
out_header = hy_obj.get_header()
out_header['bands'] = len(new_waves)
out_header['wavelength'] = new_waves.tolist()
writer = WriteENVI(spectral,out_header)
iterator =hy_obj.iterate(by = 'line')

while not iterator.complete:
    if iterator.current_line%100:
        print(iterator.current_line)
    line = iterator.read_next()[:,:(hy_obj.bands//3) * bin_dict[sensor]]
    line  = np.nanmean(view_as_blocks(line,(1,bin_dict[sensor],)),axis=(2,3))
    interpolator = interp1d(agg_waves,line,fill_value = 'extrapolate', kind = 'cubic')
    line = interpolator(new_waves)
    writer.write_line(line,iterator.current_line)












