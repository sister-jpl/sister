%matplotlib inline

import os
import hytools as ht
import pandas as pd
import isofit
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob

from isofit.core.forward import ForwardModel
from isofit.core.geometry import Geometry
from isofit.inversion.inverse import Inversion
from isofit.configs.configs import Config
from scipy.interpolate import interp1d
from isofit.utils import surface_model
import datetime as dt
import ray
from skimage.filters import threshold_otsu
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
from hytools.io.envi import WriteENVI,envi_header_dict
from scipy.signal import gaussian
import shutil

home = os.path.expanduser('~')


base_names = ['20200715103934_20200715103938_0001',
              '20210209104043_20210209104047_0001',
              '20200721104249_20200721104253_0001']
export_figures = True

for base_name in base_names:

    iso_base = '%s/isotest/isofit/' % home
    figure_dir = '%s/sister/figures/prisma/wavelength_shift/' % home
    doy = dt.datetime.strptime(base_name[:8],'%Y%m%d').timetuple().tm_yday
    base_file =  '%s/data/prisma/rad/PRISMA_%s/PRISMA_%s_' % (home, base_name,base_name)
    cal_dir = '%s/temp/%s' % (home,base_name)

    if os.path.isdir(cal_dir):
        shutil.rmtree(cal_dir)

    os.mkdir(cal_dir)
    for subdir in ['/configs','/data','/output','/lut_full/']:
        os.mkdir(cal_dir+subdir)

    #Load datasets
    radiance = ht.HyTools()
    radiance.read_file(base_file + 'rad_unprj', 'envi')
    radiance.create_bad_bands([[300,380.0],[1300.0,1500],[1780.0,1975.0],[2450.0,2600]])

    observables = ht.HyTools()
    observables.read_file(base_file + 'obs_unprj', 'envi')

    location = ht.HyTools()
    location.read_file(base_file + 'loc_unprj', 'envi')

    obs_mean  =  observables.get_chunk(0,-1,0,-1).mean(axis =(0,1))
    loc_mean  =  location.get_chunk(0,-1,0,-1).mean(axis =(0,1))

    # Export wavelength file
    wavelength_file = '%s/data/wavelengths.txt' % cal_dir
    wave_fit = pd.DataFrame(index = [x for x in range(radiance.bands)])
    wave_fit[0] =np.round(radiance.wavelengths/1000,5)
    wave_fit[1] = np.round(radiance.fwhm/1000,5)
    wave_fit.to_csv(wavelength_file, sep = '\t',header=False)

    # Create decorrelated surface model
    windows= [[380.0, 1300.0], [1500, 1780.0], [1975.0, 2450.0]]
    surface_data_path = "%s/isotest/data/surface/surface_20210215_emit" % home
    surface_file = "%s/data/surface_emit.mat" % cal_dir
    surface_config = "%s/configs/surface.json" % cal_dir
    surface_config_gen_decorr(surface_data_path,
                       windows,
                       wavelength_file,
                       surface_file,
                       surface_config)
    surface_model(surface_config)


    # Create modtran template file
    engine = 'sRTMnet'
    template_file = '%s/configs/radtran_template.json' % cal_dir
    modtran_template_create(obs_mean.astype(np.float64),
                            loc_mean.astype(np.float64),
                            doy,
                            template_file)

    # Create forward model configs
    rtm_config,surface_config,instrument_config =fw_configs(iso_base,
                                                            cal_dir,
                                                            template_file,
                                                            wavelength_file,
                                                            surface_file,
                                                            1,
                                                            "multicomponent_surface",
                                                            engine,
                                                            observables)

    full_config = Config({"forward_model":{"instrument":instrument_config,
                                           "surface":surface_config,
                                           "radiative_transfer": rtm_config}})

    if ray.is_initialized():
        ray.shutdown()
    ray.init(num_cpus = 48)

    fm = ForwardModel(full_config)

    inversion_settings = {"implementation": {"mode": "inversion",
                                             "inversion": {
                                              "windows": windows}}}
    inverse_config = Config(inversion_settings)
    iv = Inversion(inverse_config, fm)

    # Refine water vapor LUT
    water = []
    window  =5
    for line in range(10,990,100):
        print(line)
        obs_mean  =  observables.get_chunk(2,990,
                                          line-window,line+window).mean(axis =(0,1))
        loc_mean  =  location.get_chunk(2,990,
                                          line-window,line+window).mean(axis =(0,1))
        rad_mean  = radiance.get_chunk(2,990,
                                     line-window,line+window).mean(axis =(0,1))
        geom = Geometry(obs = obs_mean,loc = loc_mean)
        state_trajectory = iv.invert(rad_mean, geom)
        water.append(state_trajectory[-1][-1])

    #Create new water table
    water_min = np.min(water)*.9
    water_max = np.max(water)*1.1
    water_vals = np.linspace(water_min,water_max,5).tolist()

    # Recreate configs with refined water vapor
    for file in glob.glob(cal_dir+"/lut_full/*"):
        os.remove(file)

    rtm_config,surface_config,instrument_config =fw_configs(iso_base,
                                                            cal_dir,
                                                            template_file,
                                                            wavelength_file,
                                                            surface_file,
                                                            1,
                                                            "multicomponent_surface",
                                                            engine,
                                                            observables,
                                                            water_vals)

    full_config = Config({"forward_model":{"instrument":instrument_config,
                                           "surface":surface_config,
                                           "radiative_transfer": rtm_config}})

    fm = ForwardModel(full_config)

    inverse_config = Config(inversion_settings)
    iv = Inversion(inverse_config, fm)

    obs_mean  =  observables.get_chunk(2,990,
                                      400,600).mean(axis =(0,1))
    loc_mean  =  location.get_chunk(2,990,
                                      400,600).mean(axis =(0,1))
    rad_mean  = radiance.get_chunk(2,990,
                                      400,600).mean(axis =(0,1))

    state_trajectory = iv.invert(rad_mean, geom)

    state_est = state_trajectory[-1]
    rfl_est1 = state_est[:radiance.bands]

    rfl_smooth = np.copy(rfl_est1)
    rfl_smooth[:30]= savgol_filter(rfl_smooth[:30],7,3)
    d_start,d_end = 30,94
    rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],31,2)
    d_start,d_end = 112,140
    rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],17,4)
    d_start,d_end = 161,223
    rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],23,2)

    new_state = state_est.copy()
    new_state[:radiance.bands] =rfl_smooth
    smooth_rad = fm.calc_rdn(new_state,geom)

    rfl_est1[radiance.bad_bands] = np.nan
    rfl_smooth[radiance.bad_bands] = np.nan
    smooth_rad[radiance.bad_bands] = np.nan

    #Plot initial inversion
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    ax.plot(radiance.wavelengths,rfl_est1,c='k',lw=2)
    ax.set_xlabel("Wavelength (nm)",fontsize = 15)
    ax.set_ylabel("Reflectance",fontsize = 15)
    ax.set_xlim(375,2500)
    for direction in ['left','right','top','bottom']:
        ax.spines[direction].set_linewidth(1.5)
    if export_figures:
        plt.savefig("%s/PRISMA_%s_initial_rfl.png" % (figure_dir,base_name),
                    bbox_inches = 'tight', dpi = 500)
        plt.show()
        plt.close()

    #Plot initial inversion and smooth reflectance
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    ax.plot(radiance.wavelengths,rfl_est1,c='k',lw=2,alpha=.5)
    ax.plot(radiance.wavelengths,rfl_smooth,c='r',lw=2)
    ax.set_xlabel("Wavelength (nm)",fontsize = 15)
    ax.set_ylabel("Reflectance",fontsize = 15)
    ax.set_xlim(375,2500)
    for direction in ['left','right','top','bottom']:
        ax.spines[direction].set_linewidth(1.5)
    if export_figures:
        plt.savefig("%s/PRISMA_%s_smoothed_rfl.png" % (figure_dir,base_name),
                    bbox_inches = 'tight', dpi = 500)
        plt.show()
        plt.close()

    #Plot smoother at-sensor radiance
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    ax.plot(radiance.wavelengths,smooth_rad,c='g',lw=2)
    ax.set_xlabel("Wavelength (nm)",fontsize = 15)
    ax.set_ylabel("Radiance",fontsize = 15)
    ax.set_xlim(375,2500)
    for direction in ['left','right','top','bottom']:
        ax.spines[direction].set_linewidth(1.5)
    if export_figures:
        plt.savefig("%s/PRISMA_%s_smoothed_atsensor_radiance.png" % (figure_dir,base_name),
                    bbox_inches = 'tight', dpi = 500)
        plt.show()
        plt.close()

    shift_surf = np.zeros((radiance.lines,radiance.bands))

    #Run every 10 lines and average
    lines = [x for x in range(10,radiance.lines-10,10)]
    for line in lines:
        print(line)

        # Get mean line data
        obs_mean  =  observables.get_chunk(2,990,
                                          line-5,line+5).mean(axis =(0,1))
        loc_mean  =  location.get_chunk(2,990,
                                          line-5,line+5).mean(axis =(0,1))
        rad_mean  = radiance.get_chunk(2,990,
                                     line-5,line+5).mean(axis =(0,1))

        # Perform inversion
        geom = Geometry(obs = obs_mean,loc = loc_mean)
        state_trajectory = iv.invert(rad_mean, geom)
        state_est = state_trajectory[-1]
        rfl_est1 = state_est[:radiance.bands]
        rfl_est1[radiance.bad_bands] = np.nan

        # Generate smoothed reflectance and forward model to at-sensor radiance
        rfl_smooth = np.copy(rfl_est1)
        rfl_smooth[:30]= savgol_filter(rfl_smooth[:30],7,3)
        d_start,d_end = 30,94
        rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],31,2)
        d_start,d_end = 112,140
        rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],17,4)
        d_start,d_end = 161,223
        rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],17,2)
        d_start,d_end = 161,181
        d_start,d_end = 181,215
        new_state = state_est.copy()
        new_state[:radiance.bands] =rfl_smooth
        smooth_rad = fm.calc_rdn(new_state,geom)

        # Fit line
        new_rad = np.copy(rad_mean)
        for d_start,d_end in [[0,63],[63,94],[112,140],[161,223]]:
            mask = ~radiance.bad_bands[d_start:d_end]
            waves= radiance.wavelengths[d_start:d_end][mask]
            x = np.arange(d_end-d_start)[mask]
            smooth_rad_sub = smooth_rad[d_start:d_end]

            def poly_2(abc):
                a,b,c =abc
                shift= a + x*b + c*(x**2)
                detect_interpolator = interp1d(waves+shift,rad_mean[d_start:d_end][mask],
                                               fill_value = "extrapolate",kind='cubic')
                rad_shift = detect_interpolator(waves)
                rmse = np.sqrt(np.mean((rad_shift-smooth_rad_sub)**2))
                return rmse

            res = minimize(poly_2, (.0,.0,0), method='Nelder-Mead', tol=1e-6)
            a,b,c= res.x
            x_full = np.arange(d_end-d_start)
            shift= a + x_full*b + c*(x_full**2)
            shift_surf[line,d_start:d_end] = shift
            detect_interpolator = interp1d(waves+shift,rad_mean[d_start:d_end][mask],
                                           fill_value = "extrapolate",kind='cubic')
            new_rad[d_start:d_end]= detect_interpolator(waves)

        if line ==500:
            shift_state = iv.invert(new_rad, geom)[-1,:-2]
            shift_state[radiance.bad_bands] = np.nan
            fig = plt.figure(figsize = (6,4))
            ax = fig.add_subplot(111)
            ax.plot(radiance.wavelengths,rfl_est1,alpha=1,c = 'k', label ='Unshifted')
            ax.plot(radiance.wavelengths,shift_state,alpha=1,c = 'r',label = 'Shifted')
            ax.legend(frameon = False,loc = 'lower right',fontsize= 12)
            ax.set_xlabel("Wavelength (nm)",fontsize = 15)
            ax.set_ylabel("Reflectance",fontsize = 15)
            ax.set_xlim(375,2500)

            for direction in ['left','right','top','bottom']:
                ax.spines[direction].set_linewidth(1.5)

            if export_figures:
                plt.savefig("%s/PRISMA_%s_rfl_optimized.png" % (figure_dir,base_name),
                                bbox_inches = 'tight', dpi = 500)
                plt.show()
                plt.close()


    example_waves = [500,750,850,1200,1660,2200]
    example_bands = [radiance.wave_to_band(wave) for wave in example_waves]

    fig = plt.figure(figsize = (11,6))
    a=1
    shift_surf_smooth = np.zeros((radiance.lines,radiance.bands))
    for band in range(radiance.bands):
        wavelength= radiance.wavelengths[band]
        shift = np.copy(shift_surf[lines,band])
        z = np.polyfit(lines, shift, 3)
        p = np.poly1d(z)
        shift_surf_smooth[:,band] = p(np.arange(radiance.lines))

        if band in example_bands:
            ax = fig.add_subplot(2,3,a)
            ax.scatter(lines,shift,c='gray')
            ax.plot(shift_surf_smooth[:,band],c='r',lw=2)

            if a >3:
                ax.set_xlabel('Cross track pixel',fontsize = 15)
            if a in [1,4]:
                ax.set_ylabel('Shift (nm)',fontsize = 15)
            ax.text(.02,.9,'%s nm' % int(wavelength),transform=ax.transAxes,fontsize = 13)

            a+=1
    if export_figures:
        plt.savefig("%s/PRISMA_%s_shift_cross_track_fit.png" % (figure_dir,base_name),
                        bbox_inches = 'tight', dpi = 500)
        plt.show()
        plt.close()

    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)
    im = ax.matshow(shift_surf_smooth,
                    cmap = plt.get_cmap('seismic'),
                    vmin=-5,vmax=5)
    ax.set_aspect(.3)
    cb = plt.colorbar(im)
    cb.set_label('Shift (nm)',fontsize =15)
    ax.set_ylabel('Cross track pixel',fontsize = 15)
    if export_figures:
        plt.savefig("%s/PRISMA_%s_shift_cross_track_fit.png" % (figure_dir,base_name),
                        bbox_inches = 'tight', dpi = 500)
        plt.show()
        plt.close()



    shift_surf_smooth[:,radiance.bad_bands]=0
    shift_header = envi_header_dict()
    shift_header ['lines']= shift_surf_smooth.shape[0]
    shift_header ['samples']= shift_surf_smooth.shape[1]
    shift_header ['bands']= 1
    shift_header ['interleave']= 'bsq'
    shift_header ['data type'] = 4
    shift_header['byte order'] = 0
    shift_header['band_names'] = ['wavelength shift(nm)']

    shift_file = "%s/sister/data/prisma/wavelength_shift/PRISMA_%s_wavelength_shift" % (home,base_name)
    writer = WriteENVI(shift_file,shift_header )
    writer.write_band(shift_surf_smooth,0)












