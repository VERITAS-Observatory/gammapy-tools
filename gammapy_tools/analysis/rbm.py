import logging
import yaml
import numpy as np
from scipy.stats import norm
import os
from astropy.io import fits
from os import environ
from astropy.table import Table
# from scipy.optimize import fsolve
from typing import Optional, Tuple
from tqdm.auto import tqdm

# %matplotlib inline
import astropy.units as u
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion
import matplotlib.pyplot as plt

import gammapy
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import MapDatasetOnOff,MapDataset
from gammapy.estimators import ExcessMapEstimator
from gammapy.makers import RingBackgroundMaker,SafeMaskMaker,MapDatasetMaker
from gammapy.maps.wcs.ndmap import WcsNDMap
from gammapy.modeling.models import SpectralModel
from gammapy.data import DataStore
from gammapy.maps import MapAxis, WcsGeom

import seaborn as sns
sns.set_theme(font="MathJax_Main",font_scale=2,style='ticks',context='paper',palette='pastel')

# importing package version
from .._version import __version__ as gpt_version

__notebook_version__ = "0.3"

log = logging.getLogger(__name__)


# No longer needed:
#
# def estimate_alpha(S: float, N_on: float, N_off: float) -> float:
#     """
#     Numerically estimates alpha from significance, ON counts, and OFF counts
# 
#     Parameters
#     ----------
#         S: significance
#         N_on: on counts
#         N_off: off counts*alpha
# 
#     Returns
#     ----------
#         alpha: estimated alpha
#     """
# 
#     # Define the function to find the root of
#     def equation(alpha):
#         # This is Li & Ma equation 17
#         return (
#             np.sqrt(2)
#             * np.sqrt(
#                 N_on * np.log((1 + alpha) / alpha * (N_on / (N_on + (N_off / alpha))))
#                 + (N_off / alpha)
#                 * np.log((1 + alpha) * ((N_off / alpha) / (N_on + (N_off / alpha))))
#             )
#             - S
#         )
# 
#     # Initial guess for alpha
#    alpha_initial_guess = 1e-2
# 
#     # Solve for alpha
#     alpha_solution = fsolve(equation, alpha_initial_guess)
#     return alpha_solution[0]


def rbm_analysis(
    config: dict,
) -> Tuple[float, float, float, float, np.ndarray, float, np.ndarray, np.ndarray, np.ndarray]:
    """
    Performs a basic RBM analysis

    Parameters
    ----------
        config: configuration file

    Returns
    ----------
        counts: total counts
        background: off counts
        alpha: relative size of the on/off regions & exposure times
        sigma: significance at the source location (defined in the config file)
        excess_map: map of excess counts
        exposure: time on source
        significance_map: significance map
        exclusion_mask: exclusion mask
        alpha_map: alpha map
    """
    if not os.path.exists(config["io"]["results_dir"]):
        os.makedirs(config["io"]["results_dir"])

    data_store = DataStore.from_dir(config["io"]["out_dir"])
    
    # select only observations from runlist, if specified
    if config["io"]["from_runlist"]:
        obs_ids = np.genfromtxt(
            config["io"]["runlist"], unpack=True
        ).tolist()
        observations = data_store.get_observations(obs_id=obs_ids,required_irf="full-enclosure")
    else:
        observations = data_store.get_observations(required_irf="full-enclosure")

    if config["run_selection"]["pos_from_DL3"]:
        # get RA and DEC from first run
        hdul = fits.open(
            config["io"]["out_dir"] + os.listdir(config["io"]["out_dir"])[0]
        )
        source_pos = SkyCoord(
            hdul[1].header["RA_OBJ"] * u.deg, hdul[1].header["DEC_OBJ"] * u.deg
        )

    else:
        source_pos = SkyCoord(
            config["run_selection"]["source_ra"],
            config["run_selection"]["source_dec"],
            frame="icrs",
            unit="deg",
        )

    energy_axis = MapAxis.from_energy_bounds(
        config["sky_map"]["e_min"],config["sky_map"]["e_max"],10,unit="TeV"
    )
    # Reduced IRFs are defined in true energy (i.e. not measured energy).
    energy_axis_true = MapAxis.from_energy_bounds(
        config["sky_map"]["e_min"], config["sky_map"]["e_max"], 200, unit="TeV", name="energy_true"
    )

    map_deg = config["sky_map"]["map_deg"]
    geom = WcsGeom.create(
        skydir=(source_pos.ra.value, source_pos.dec.value),
        binsz=config["sky_map"]["bin_sz"],
        width=(map_deg, map_deg),
        frame="fk5",
        proj="TAN",
        axes=[energy_axis],
    )
    geom_image = geom.to_image().to_cube([energy_axis.squash()])

    # Make the exclusion mask
    on_region = CircleSkyRegion(
        center=source_pos, radius=config["sky_map"]["on_exclusion_region"] * u.deg
    )
    all_ex = [on_region]
    if len(config["sky_map"]["exclusion_regions"]) != 0:
        for region in config["sky_map"]["exclusion_regions"]:
            ra, dec = region[0]
            radius = region[1]
            all_ex.append(
                CircleSkyRegion(
                    center=SkyCoord(ra, dec, unit="deg", frame="icrs"),
                    radius=radius * u.deg,
                )
            )

    star_data = np.loadtxt(
        # environ["GAMMAPY_DATA"] + "/catalogs/Hipparcos_MAG8_1997.dat", usecols=(0, 1, 2, 3, 4)
        environ["GAMMAPY_DATA"] + "/catalogs/Hipparcos_MAG8_1997.dat",
        usecols=(0, 1, 2, 3),
    )
    star_cat = Table(
        {
            "ra": star_data[:, 0],
            "dec": star_data[:, 1],
            "id": star_data[:, 2],
            "mag": star_data[:, 3],
            # "colour": star_data[:, 4],
        }
    )
    star_mask = (
        np.sqrt(
            (star_cat["ra"] - source_pos.ra.deg) ** 2
            + (star_cat["dec"] - source_pos.dec.deg) ** 2
        )
        < 2.0
    )
    star_mask &= (star_cat["mag"]) < 8

    for src in star_cat[star_mask]:
        all_ex.append(
            CircleSkyRegion(
                center=SkyCoord(src["ra"], src["dec"], unit="deg", frame="icrs"),
                radius=0.3 * u.deg,
            )
        )

    exclusion_mask = ~geom_image.region_mask(all_ex)
    # exclusion_mask.sum_over_axes().plot()
    
    stacked = MapDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name="stacked"
    )

    # the only safemaskmaker parameter we want for sky maps is offset-max
    offset_max = config["sky_map"]["offset_max"] * u.deg
    maker = MapDatasetMaker()
    maker_safe_mask = SafeMaskMaker(
        methods=["offset-max"], offset_max=offset_max
    )

    # create stacked dataset with all observations
    for obs in tqdm(observations):
        dataset = maker.run(stacked,obs)
        dataset = maker_safe_mask.run(dataset, obs)
        stacked.stack(dataset)

    ring_maker = RingBackgroundMaker(
        r_in=config["sky_map"]["ring_rin"],
        width=config["sky_map"]["ring_width"],
        exclusion_mask=exclusion_mask,
    )
    dataset_on_off = ring_maker.run(stacked)
    
    output_dataset = dataset_on_off.to_spectrum_dataset(
        CircleSkyRegion(center=source_pos, radius=config["sky_map"]["theta"] * u.deg),
        containment_correction=True,
    )
    output_dict = output_dataset.info_dict()

    estimator = ExcessMapEstimator(
        config["sky_map"]["theta"] * u.deg,
        selection_optional=["alpha"],
        # spectral_model=spectral_model,
        correlate_off=False,
    )
    lima_maps = estimator.run(dataset_on_off)
    significance_map = lima_maps["sqrt_ts"]
    excess_map = lima_maps["npred_excess"]
    alpha_map = lima_maps["alpha"]

    counts = lima_maps["npred"].get_by_coord(
        [source_pos.ra, source_pos.dec, 1 * u.TeV]
    )[0]
    background = lima_maps["npred_background"].get_by_coord(
        [source_pos.ra, source_pos.dec, 1 * u.TeV]
    )[0]
    sigma = lima_maps["sqrt_ts"].get_by_coord(
        [source_pos.ra, source_pos.dec, 1 * u.TeV]
    )[0]
    alpha = alpha_map.get_by_coord(
        [source_pos.ra, source_pos.dec, 1 * u.TeV]
    )[0]
    exposure = output_dict["ontime"]
    
    return (
        counts,
        background,
        alpha,
        sigma,
        excess_map,
        exposure,
        significance_map,
        exclusion_mask,
        alpha_map,
    )


def rbm_plots(
    config: dict,
    spectral_points: Table,
    excess_map: WcsNDMap,
    significance_map: WcsNDMap,
    exclusion_mask: WcsNDMap,
    alpha_map: WcsNDMap,
    save: Optional[bool] = True,
    plot: Optional[bool] = True,
    spectrum: Optional[bool] = True,
) -> None:
    """
    Makes + optionally saves significance/excess maps,
    significance distribution, and cumulative significance
    Parameters
    ----------
        config: configuration file
        spectral_points: spectral points
        excess map: excess counts map
        significance map: Li&Ma significance map
        exclusion mask: exclusion mask of regions to be excluded
        alpha_map: map of Li&Ma alpha values
        save: if True, saves plots with prefix defined in config['plot_names']
    """

    # significance & excess plots
    fig, (ax1, ax2, ax3) = plt.subplots(
        figsize=(18, 5), subplot_kw={"projection": significance_map.geom.wcs}, ncols=3
    )
    cmap=sns.color_palette("magma", as_cmap=True)
    ax1.set_title("Significance map")
    if config["sky_map"]["truncate"]:
        significance_map.plot(ax=ax1, add_cbar=True, vmin=-5, vmax=5,cmap=cmap)
    else:
        significance_map.plot(ax=ax1, add_cbar=True,cmap=cmap)

    ax2.set_title("Excess map")
    excess_map.plot(ax=ax2, add_cbar=True,cmap=cmap)

    ax3.set_title("Alpha map")
    alpha_map.plot(ax=ax3, add_cbar=True, cmap=cmap)
    
    plt.savefig(
        config["plot_names"] + "sig_excess_alpha.png", format="png", bbox_inches="tight"
    )
    plt.show()

    # create a 2D mask for the images
    # significance_map_off = significance_map * exclusion_mask
    significance_all = significance_map.data[
        np.isfinite(significance_map.data)
    ].flatten()
    significance_off = significance_map.data[
        exclusion_mask & np.isfinite(significance_map.data)
    ].flatten()

    fig, ax = plt.subplots()
    ax.hist(
        significance_all,
        density=True,
        alpha=0.9,
        color="palevioletred",
        label="all bins",
        bins=np.linspace(-5, 10, 50),
    )

    ax.hist(
        significance_off,
        density=True,
        alpha=0.9,
        color="darkseagreen",
        label="off bins",
        bins=np.linspace(-5, 10, 50),
    )

    # Now, fit the off distribution with a Gaussian
    mu, std = norm.fit(significance_off)
    x = np.linspace(-10, 10, 100)
    p = norm.pdf(x, mu, std)
    gauss = norm.pdf(x,0,1)
    ax.plot(x,gauss,lw=2,ls='--',color='black')
    ax.plot(x, p, lw=2, color="black")
    ax.legend()
    ax.set_xlabel("Significance")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    # unused values
    # xmin, xmax = np.min(significance_all), np.max(significance_all)
    ax.set_xlim(-5, 10)

    print(f"Fit results: mu = {mu:.2f}, std = {std:.2f}")
    ax.text(-4.5, 0.5, f"Fit results: mu = {mu:.2f}, std = {std:.2f}")
    if plot:
        plt.show()
    else:
        plt.clf()

    if spectrum:
        spectral_points.plot(sed_type="dnde")
        plt.savefig(
            config["io"]["results_dir"] + config["plot_names"] + "spectral_points.png",
            format="png",
            bbox_inches="tight",
        )
        if plot:
            plt.show()

    return


def write_validation_info(
    config: dict,
    spectral_model: SpectralModel,
    counts: float,
    background: float,
    alpha: float,
    sigma: float,
    exposure: float,
    spectrum: Optional[bool] = True
) -> None:

    if not os.path.exists(config["io"]["results_dir"]):
        os.makedirs(config["io"]["results_dir"])

    if spectrum:
        # calculate integral flux from spectral model
        e_min = config["spectrum"]["e_min"]*u.TeV
        e_max = config["spectrum"]["e_max"]*u.TeV
        
        flux,flux_err = spectral_model.integral_error(e_min,e_max)
        spectab = spectral_model.to_dict()
        if 'components' in spectab.keys():
            spectab = spectab['components'][0]
        index = spectab["spectral"]["parameters"][0]["value"]
        index_err = spectab["spectral"]["parameters"][0]["error"]
        norm = spectab["spectral"]["parameters"][1]["value"]
        norm_err = spectab["spectral"]["parameters"][1]["error"]
    
    else:
        flux=flux_err=index=index_err=norm=norm_err = 0
    output_dict = {
        "analysis notebook version": __notebook_version__,
        "gammapy-tools version": gpt_version,
        "source": config["run_selection"]["source_name"],
        "gammapy version": gammapy.__version__,
        "exposure (min)": float(exposure.value) / 60,
        "on": int(counts),
        "off": f"{background:.2f}",
        "alpha": f"{alpha:.3e}",
        "significance": f"{sigma:.2f}",
        "integral_flux": f"{flux:.2e}",
        "flux_err": f"{flux_err:.2e}",
        "flux_UL": "n/a",
        "norm": f"{norm:.2e}",
        "norm_err": f"{norm_err:.2e}",
        "index": f"{index:.2f}",
        "index_err": f"{index_err:.2f}",
    }

    with open(config["io"]["results_dir"] + config["results_file"], "w") as outfile:
        yaml.dump(output_dict, outfile, sort_keys=False)

    print("======== RESULTS ========")
    for key, value in output_dict.items():
        print(key, ":", value)
    
    # print in wiki format
    print(f'\n||notebook v2.0||{float(exposure.value) / 60}||{int(counts)}||{int(background)}||{alpha:.3e}||{sigma:.2f}||{flux:.2e} cm2 s-1||{flux_err:.2e} cm2 s-1||{norm:.2e} TeV-1 s-1 cm-2||{norm_err:.2e} TeV-1 s-1 cm-2||{index:.2f}||{index_err:.2f}||')
    return
