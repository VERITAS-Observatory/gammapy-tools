import logging
import yaml
import numpy as np
from scipy.stats import norm
from IPython.display import display
import os
from astropy.io import fits

# %matplotlib inline
import astropy.units as u
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion
import matplotlib.pyplot as plt

import gammapy
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import MapDatasetOnOff
from gammapy.estimators import ExcessMapEstimator
from gammapy.makers import RingBackgroundMaker
from gammapy.data import DataStore

from gammapy.modeling.models import PowerLawSpectralModel

from astroquery.simbad import Simbad

log = logging.getLogger(__name__)


def rbm_analysis(config):
    """
    Performs a basic RBM analysis

    Parameters
    ----------
        config: configuration file

    Returns
    ----------
        counts: total counts
        excess: excess counts
        background: off counts
        alpha: relative size of the on/off regions & exposure times
        sigma: significance at the source location (defined in the config file)
        excess_map: map of excess counts
        significance_map: significance map
        significance_map_off: background significance map
    """

    data_store = config["io"]["out_dir"]

    map_deg = config["sky_map"]["map_deg"]
    source_config = AnalysisConfig()
    source_config.datasets.type = "3d"
    source_config.observations.datastore = data_store
    
    #select only observations from runlist, if specified
    if config["io"]["from_runlist"]:
        source_config.observations.obs_ids = np.genfromtxt(config["io"]["runlist"],unpack=True).tolist()
    
    if config["run_selection"]["pos_from_DL3"]:
        #get RA and DEC from first run
        hdul = fits.open(config["io"]["out_dir"]+os.listdir(config["io"]["out_dir"])[0])
        source_pos = SkyCoord(hdul[1].header["RA_OBJ"]*u.deg, hdul[1].header["DEC_OBJ"]*u.deg)
    
    else:
        source_pos = SkyCoord(
            config["run_selection"]["source_ra"],
            config["run_selection"]["source_dec"],
            frame="icrs",
            unit="deg",
        )
    
    source_config.datasets.geom.wcs.skydir = {
        "lon": source_pos.ra,
        "lat": source_pos.dec,
        "frame": "icrs",
    }

    source_config.observations.required_irf = ["aeff", "edisp"]

    source_config.datasets.geom.wcs.width = {
        "width": f"{map_deg} deg",
        "height": f"{map_deg} deg",
    }
    source_config.datasets.geom.wcs.binsize = config["sky_map"]["bin_size"]
    source_config.datasets.map_selection = ["counts", "exposure", "background", "edisp"]

    # Cutout size (for the run-wise event selection)
    source_config.datasets.geom.selection.offset_max = map_deg * u.deg

    # We now fix the energy axis for the counts map - (the reconstructed energy binning)
    source_config.datasets.geom.axes.energy.min = "0.1 TeV"
    source_config.datasets.geom.axes.energy.max = "100 TeV"
    source_config.datasets.geom.axes.energy.nbins = 30

    source_config.excess_map.correlation_radius = str(config["sky_map"]["theta"]) +  " deg"
    
    # We need to extract the ring for each observation separately, hence, no stacking at this stage
    source_config.datasets.stack = False

    source_config.datasets.safe_mask.parameters = {'aeff_percent':config["sky_map"]["aeff_max_percent"], 'offset_max':config["sky_map"]["offset_max"]*u.deg}
    source_config.datasets.safe_mask.methods = ['aeff-max','offset-max']

    analysis = Analysis(source_config)

    analysis.config.datasets.geom.axes.energy_true = (
        analysis.config.datasets.geom.axes.energy
    )

    analysis.get_observations()
    analysis.get_datasets()

    simbad = Simbad()
    simbad.reset_votable_fields()
    simbad.add_votable_fields("ra", "dec", "flux(B)", "flux(V)", "jp11")
    simbad.remove_votable_fields("coordinates")

    srcs_tab = simbad.query_region(source_pos, radius=1.5 * u.deg)
    srcs_tab = srcs_tab[srcs_tab["FLUX_B"] < config["sky_map"]["min_star_brightness"]]
    srcs_tab = srcs_tab[srcs_tab["FLUX_V"] != np.ma.masked]

    # get the geom that we use
    geom = analysis.datasets[0].counts.geom
    energy_axis = analysis.datasets[0].counts.geom.axes["energy"]
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
    stars = []
    for star in srcs_tab:
        pos = SkyCoord(star["RA"], star["DEC"], frame="fk5", unit=(u.hourangle, u.deg))
        star = CircleSkyRegion(center=pos, radius=0.3 * u.deg)
        stars.append(star)
        all_ex.append(star)
    exclusion_mask = ~geom_image.region_mask(all_ex)
    # exclusion_mask.sum_over_axes().plot()

    ring_maker = RingBackgroundMaker(
        r_in=config["sky_map"]["ring_rin"],
        width=config["sky_map"]["ring_width"],
        exclusion_mask=exclusion_mask,
    )

    energy_axis_true = analysis.datasets[0].exposure.geom.axes["energy_true"]
    stacked_on_off = MapDatasetOnOff.create(
        geom=geom_image, energy_axis_true=energy_axis_true, name="stacked"
    )

    for dataset in analysis.datasets:
        # Ring extracting makes sense only for 2D analysis
        dataset_on_off = ring_maker.run(dataset.to_image())
        stacked_on_off.stack(dataset_on_off)

    # we'll convert our map dataset into a spectrum dataset,
    # which is more helpful for extracting info
    output_dataset = stacked_on_off.to_spectrum_dataset(
        CircleSkyRegion(center=source_pos, radius=config["sky_map"]["theta"] * u.deg),
        containment_correction=True,
    )
    output_dict = output_dataset.info_dict()

    # spectral model for estimator
    amp, idx = config["spectrum"]["params"]
    spectral_model = PowerLawSpectralModel(
        amplitude=float(amp) * u.Unit("cm-2 s-1 TeV-1"),
        index=float(idx),
        reference=1 * u.TeV,
    )

    estimator = ExcessMapEstimator(
        config["sky_map"]["theta"] * u.deg,
        selection_optional=["all"],
        spectral_model=spectral_model,
    )
    lima_maps = estimator.run(stacked_on_off)
    significance_map = lima_maps["sqrt_ts"]
    excess_map = lima_maps["npred_excess"]

    counts = output_dict["counts"]
    excess = output_dict["excess"]
    background = output_dict["background"]
    alpha = output_dict["alpha"]
    sigma = output_dict["sqrt_ts"]
    exposure = output_dict["ontime"]

    #significance_map_off = significance_map * exclusion_mask
    #significance_map_off = significance_map[exclusion_mask]

    return (
        counts,
        excess,
        background,
        alpha,
        sigma,
        excess_map,
        exposure,
        significance_map,
        significance_map_off,
        exclusion_mask,
    )


def rbm_plots(
    config,
    spectral_points,
    excess_map,
    significance_map,
    exclusion_mask,
    c_sig,
    c_time,
    save=True,
    plot=True,
):
    """
    Makes + optionally saves significance/excess maps,
    significance distribution, and cumulative significance
    Parameters
    ----------
        config: configuration file
        spectral_points: spectral points
        excess map: excess counts map
        significance map: Li&Ma significance map
        exclusion mask: Exclusion mask of regions to be excluded
        c_sig: cumulative significance
        c_time: cumulative time
        save: if true, saves plots with prefix defined in config['plot_names']
    """

    # significance & excess plots
    fig, (ax1, ax2) = plt.subplots(
        figsize=(11, 5), subplot_kw={"projection": significance_map.geom.wcs}, ncols=2
    )

    ax1.set_title("Significance map")
    if config["sky_map"]["truncate"]:
        significance_map.plot(ax=ax1, add_cbar=True, vmin=-5, vmax=5)
    else:
        significance_map.plot(ax=ax1, add_cbar=True)

    ax2.set_title("Excess map")
    excess_map.plot(ax=ax2, add_cbar=True)
    plt.savefig(
        config["plot_names"] + "sig_excess.png", format="png", bbox_inches="tight"
    )
    plt.show()

    # significance distribution
    significance_all = significance_map.data[np.isfinite(significance_map.data)]
    significance_off = significance_map.data[
        np.isfinite(significance_map.data) & exclusion_mask
    ]

    fig, ax = plt.subplots()
    ax.hist(
        significance_all,
        density=True,
        alpha=0.5,
        color="red",
        label="all bins",
        bins=np.linspace(-5, 10, 100),
    )

    ax.hist(
        significance_off,
        density=True,
        alpha=0.5,
        color="blue",
        label="off bins",
        bins=np.linspace(-5, 10, 100),
    )

    # Now, fit the off distribution with a Gaussian
    mu, std = norm.fit(significance_off)
    x = np.linspace(-8, 8, 50)
    p = norm.pdf(x, mu, std)
    ax.plot(x, p, lw=2, color="black")
    ax.legend()
    ax.set_xlabel("Significance")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    # xmin, xmax = np.min(significance_all), np.max(significance_all)
    ax.set_xlim(-5, 10)

    print(f"Fit results: mu = {mu:.2f}, std = {std:.2f}")
    ax.text(-4.5, 0.5, f"Fit results: mu = {mu:.2f}, std = {std:.2f}")
    plt.savefig(
        config["io"]["results_dir"] + config["plot_names"] + "sig_dist.png",
        format="png",
        bbox_inches="tight",
    )
    if plot:
        plt.show()
    else:
        plt.clf()
    # cumulative significance
    plt.plot(c_time, c_sig, "ko")
    plt.xlabel("Time [h]")
    plt.ylabel(r"Significance [$\sigma$]")
    plt.title("Cumulative Significance")
    plt.savefig(
        config["io"]["results_dir"] + config["plot_names"] + "cumulative_sig.png",
        format="png",
        bbox_inches="tight",
    )
    if plot:
        plt.show()
    else:
        plt.clf()

    # spectral points
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
    config, spectral_model, flux, counts, background, alpha, sigma, exposure
):
    fluxtab = flux.to_table(sed_type="flux", format="lightcurve")
    f200 = fluxtab["flux"].value[0][0]
    f200_err = fluxtab["flux_err"].value[0][0]
    ul = fluxtab["is_ul"].value[0][0]

    spectab = spectral_model.to_parameters_table()
    index = spectab["value"][0]
    index_err = spectab["error"][0]
    norm = spectab["value"][1]
    norm_err = spectab["error"][1]

    output_dict = {
        "analysis notebook version": 1.0,
        "source": config["run_selection"]["source_name"],
        "gammapy version": gammapy.__version__,
        "exposure": float(exposure.value) / 60,
        "on": int(counts),
        "off": int(background),
        "alpha": f"{alpha:.3e}",
        "significance": f"{sigma:.2f}",
        "flux": f"{f200:.2e}",
        "flux_err": f"{f200_err:.2e}",
        "flux_UL": bool(ul),
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

    return
