import argparse
import os
import subprocess
from isofit.utils import surface_model
from sister.utils.isofit import gen_wavelength_file,surface_config_gen,get_surface_spectra
from sister.utils.misc import download_file
import yaml

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('rdn_file', type=str,
                        help='Radiance')
    parser.add_argument('obs_file', type=str,
                        help='Observables')
    parser.add_argument('loc_file', type=str,
                         help='Location')
    parser.add_argument('output_dir', type=str,
                         help='Output directory')
    parser.add_argument('temp_dir', type=str,
                         help='Temp directory')
    parser.add_argument('config', type=str,
                        help='ISOFIT config files')
    args = parser.parse_args()

    base_name = os.path.basename(args.rdn_file)[:38]
    date = base_name[4:12]

    # Load ISOFIT config file
    with open(args.config,'r') as file:
        configs = yaml.load(file, Loader=yaml.FullLoader)

    for key in configs['surface']:
        get_surface_spectra(args.temp_dir + key,
                            configs['surface'][key])


    surface_config = "%s/surface_config.json" % args.temp_dir
    surface_file = '%s/surface_filtered.mat'% args.temp_dir
    wavelength_file = gen_wavelength_file(args.rdn_file,
                                           args.temp_dir)

    surface_config_gen(args.temp_dir,
                       configs['windows'],
                       wavelength_file,
                       surface_file,
                       surface_config)

    surface_model(surface_config)

    apply_oe  = ['python']
    apply_oe.append('%s/isofit/utils/apply_oe.py' % configs['base'])
    apply_oe.append(args.rdn_file)
    apply_oe.append(args.loc_file)
    apply_oe.append(args.obs_file)
    apply_oe.append(args.output_dir)
    apply_oe.append('NA-%s' % date)
    apply_oe += ['--surface_path', surface_file]
    apply_oe += ['--n_cores',str(configs['n_cores'])]
    apply_oe += ['--empirical_line',str(configs['empirical_line'])]
    apply_oe += ['--presolve',str(configs['presolve'])]
    apply_oe += ['--ray_temp_dir',args.temp_dir]
    apply_oe += ['--log_file','%s/%s_logfile' % (args.temp_dir,base_name)]
    apply_oe += ['--emulator_base',configs['emulator']]

    if configs['radiance_factors']:
        rad_factors = "%s/radiance_factors.txt" % args.temp_dir
        download_file(rad_factors,configs['radiance_factors'])
        apply_oe += ['--rdn_factors_path',rad_factors]

    apply = subprocess.Popen(apply_oe)
    apply.wait()

if __name__ == "__main__":
    main()