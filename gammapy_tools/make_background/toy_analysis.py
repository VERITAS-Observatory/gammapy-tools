import logging
import os

# %matplotlib inline
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import MapDatasetOnOff
from gammapy.estimators import ExcessMapEstimator
from gammapy.makers import RingBackgroundMaker
from regions import CircleSkyRegion
from scipy.stats import norm

log = logging.getLogger(__name__)


def run_rbm(driver_file):
    print(driver_file["io"]["out_dir"])
    print(driver_file["run_selection"]["source_name"])
    # source_pos = SkyCoord.from_name(driver_file["run_selection"]["source_name"])
    source_pos = SkyCoord(
        driver_file["run_selection"]["source_ra"],
        driver_file["run_selection"]["source_dec"],
        unit="deg",
    )
    # source_pos = SkyCoord(228.32, -59.08, unit="deg")

    config = AnalysisConfig()
    # Select observations - 2.5 degrees from the source position
    config.observations.datastore = driver_file["io"]["out_dir"]
    config.observations.obs_cone = {
        "frame": "icrs",
        "lon": source_pos.ra,
        "lat": source_pos.dec,
        "radius": 2.5 * u.deg,
    }

    config.datasets.type = "3d"
    config.datasets.geom.wcs.skydir = {
        "lon": source_pos.ra,
        "lat": source_pos.dec,
        "frame": "icrs",
    }  # The WCS geometry - centered on MSH 15-52

    config.datasets.geom.wcs.width = {"width": "3.5 deg", "height": "3.5 deg"}
    config.datasets.geom.wcs.binsize = "0.02 deg"

    config.datasets.safe_mask.parameters = {
        "aeff_percent": 0.15,
        "offset_max": 1.7 * u.deg,
    }
    config.datasets.safe_mask.methods = ["aeff-max", "offset-max"]

    # Cutout size (for the run-wise event selection)
    # config.datasets.geom.selection.offset_max = driver_file["binning"]["OffMax"] * u.deg

    # We now fix the energy axis for the counts map - (the reconstructed energy binning)
    config.datasets.geom.axes.energy.min = "0.1 TeV"
    config.datasets.geom.axes.energy.max = "300 TeV"
    config.datasets.geom.axes.energy.nbins = 30

    # We need to extract the ring for each observation separately, hence, no stacking at this stage
    config.datasets.stack = False

    analysis = Analysis(config)

    # for this specific case,w e do not need fine bins in true energy
    analysis.config.datasets.geom.axes.energy_true = (
        analysis.config.datasets.geom.axes.energy
    )

    # `First get the required observations
    analysis.get_observations()

    analysis.get_datasets()

    # get the geom that we use
    geom = analysis.datasets[0].counts.geom
    energy_axis = analysis.datasets[0].counts.geom.axes["energy"]
    geom_image = geom.to_image().to_cube([energy_axis.squash()])

    # Make the exclusion mask

    this_dir, this_filename = os.path.split(__file__)
    try:
        star_path = os.path.join(this_dir, "Hipparcos_MAG8_1997.dat")
        star_data = np.loadtxt(star_path, usecols=(0, 1, 2, 3), skiprows=62)
    except Exception:
        star_path = os.path.join(
            os.environ.get("GAMMAPY_DATA"), "catalogs/", "Hipparcos_MAG8_1997.dat"
        )
        star_data = np.loadtxt(star_path, usecols=(0, 1, 2, 3), skiprows=62)

    star_cat = Table(
        {
            "ra": star_data[:, 0],
            "dec": star_data[:, 1],
            "id": star_data[:, 2],
            "mag": star_data[:, 3],
        }
    )
    mask = (
        np.sqrt(
            (star_cat["ra"] - source_pos.ra.deg) ** 2
            + (star_cat["dec"] - source_pos.dec.deg) ** 2
        )
        < 2
    )
    mask &= star_cat["mag"] < 7
    srcs_tab = star_cat[mask]

    regions = CircleSkyRegion(center=source_pos, radius=0.35 * u.deg)
    excluded_regions = [regions]
    stars = []
    for star in srcs_tab:
        pos = SkyCoord(star["ra"], star["dec"], frame="fk5", unit=(u.deg, u.deg))
        star = CircleSkyRegion(center=pos, radius=0.35 * u.deg)
        stars.append(star)
        excluded_regions.append(star)

    exclusion_mask = ~geom_image.region_mask(excluded_regions)
    exclusion_mask.sum_over_axes().plot()

    ring_maker = RingBackgroundMaker(
        r_in="0.6 deg", width="0.2 deg", exclusion_mask=exclusion_mask
    )

    energy_axis_true = analysis.datasets[0].exposure.geom.axes["energy_true"]
    stacked_on_off = MapDatasetOnOff.create(
        geom=geom_image, energy_axis_true=energy_axis_true, name="stacked"
    )

    for dataset in analysis.datasets:
        # Ring extracting makes sense only for 2D analysis
        dataset_on_off = ring_maker.run(dataset.to_image())
        stacked_on_off.stack(dataset_on_off)

    # Using a convolution radius of 0.04 degrees
    estimator = ExcessMapEstimator(
        np.sqrt(driver_file["analysis_selection"]["theta2"]) * u.deg,
        selection_optional=[],
    )
    lima_maps = estimator.run(stacked_on_off)

    significance_map = lima_maps["sqrt_ts"]
    excess_map = lima_maps["npred_excess"]

    # We can plot the excess and significance maps
    fig, (ax1, ax2) = plt.subplots(
        figsize=(11, 5), subplot_kw={"projection": lima_maps.geom.wcs}, ncols=2
    )

    ax1.set_title("Significance map")
    significance_map.plot(ax=ax1, add_cbar=True)

    ax2.set_title("Excess map")
    excess_map.plot(ax=ax2, add_cbar=True)
    fig.savefig(driver_file["io"]["results_dir"] + "/skymaps.png")

    # create a 2D mask for the images
    significance_map_off = significance_map * exclusion_mask
    significance_all = significance_map.data[np.isfinite(significance_map.data)]
    significance_off = significance_map_off.data[np.isfinite(significance_map_off.data)]

    fig, ax = plt.subplots()
    binning = np.linspace(-5, 10)
    ax.hist(
        significance_all,
        density=True,
        alpha=0.5,
        color="red",
        label="all bins",
        bins=binning,
    )

    ax.hist(
        significance_off,
        density=True,
        alpha=0.5,
        color="blue",
        label="off bins",
        bins=binning,
    )

    # Now, fit the off distribution with a Gaussian
    mu, std = norm.fit(significance_off)
    x = np.linspace(-5, 5, 50)
    p = norm.pdf(x, mu, std)
    p_norm = norm.pdf(x, 0, 1)
    ax.plot(
        x, p, lw=2, color="black", label=f"Fit results: mu = {mu:.2f}, std = {std:.2f}"
    )
    ax.plot(
        x,
        p_norm,
        lw=2,
        ls="--",
        color="cyan",
        label=f"Fit results: mu = {0}, std = {1}",
    )
    ax.legend()
    ax.set_xlabel("Significance")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    ax.grid()
    ax.legend()
    # print(f"Fit results: mu = {mu:.2f}, std = {std:.2f}")
    fig.savefig(driver_file["io"]["results_dir"] + "/sigdist.png")
