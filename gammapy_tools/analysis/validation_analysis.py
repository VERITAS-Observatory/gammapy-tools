# import sys

# import yaml
# from matplotlib import pyplot as plt
# from tqdm.auto import tqdm

from .data_products import make_spectrum_RE, get_flux_lc
from .rbm import rbm_analysis, rbm_plots, write_validation_info

# sys.path.insert(1,'../background_from_bkg/makeBackground/')
# from makeBackground import *


# import warnings
# warnings.filterwarnings("ignore")

# config_path = sys.argv[1] #path to config file

# with open(config_path, 'r') as file:
#     config = yaml.safe_load(file)


def validation_analysis(config):

    # Never used...
    # data_store = DataStore.from_dir(config["io"]["in_dir"])
    # if config['io']['from_runlist']:
    #     obs_ids = np.genfromtxt(config['io']['runlist'],dtype=int,unpack=True)

    # else:
    #     obs_ids = data_store.obs_ids

    # spectral analysis
    spectral_points, spectral_model, cumulative_time, cumulative_sig = make_spectrum_RE(
        config
    )

    # intregral flux calculations
    # to get a LC, add type='runwise' for runwise or type='custom' for bin size set in config file
    flux = get_flux_lc(config, "flux")

    # RBM analysis
    (
        counts,
        excess,
        background,
        alpha,
        sigma,
        excess_map,
        exposure,
        significance_map,
        significance_map_off,
    ) = rbm_analysis(config)

    # generate plots
    rbm_plots(
        config,
        excess_map,
        significance_map,
        significance_map_off,
        cumulative_sig,
        cumulative_time,
    )

    # write analysis output for validation
    write_validation_info(
        config, spectral_model, flux, counts, background, alpha, sigma, exposure
    )
