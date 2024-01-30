import numpy as np
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion
from astropy.time import Time
from astropy.table import Table

# %matplotlib inline
import matplotlib.pyplot as plt

# gammapy imports
from gammapy.data import DataStore
from gammapy.datasets import (
    Datasets,
    FluxPointsDataset,
    SpectrumDataset,
)
from gammapy.estimators import FluxPointsEstimator
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    PowerLawSpectralModel,
    LogParabolaSpectralModel,
    SkyModel,
)
from gammapy.estimators import LightCurveEstimator


def make_spectrum_RE(config, plot=True):
    """Make a RE spectrum

    Parameters
    ----------
    observations: list of gammapy Observations object
    source_pos: astropy coordinates object containing source coordinates
    config: config file

    Returns
    ----------
    flux_points: spectral flux points object
    spectral_model: best-fit spectral model object
    time: time for cumulative significance
    sig: sqrt(ts) for cumulative significance
    """

    Emin = config["spectrum"]["Emin"]
    Emax = config["spectrum"]["Emax"]
    spectrum = config["spectrum"]["type"]
    spectral_params = config["spectrum"]["params"]
    nbins = config["spectrum"]["nbins"]
    theta = config["sky_map"]["theta"]

    datastore = DataStore.from_dir(config["io"]["out_dir"])
    observations = datastore.get_observations()

    if config["source"]["use_name"]:  # get position from name
        source_pos = SkyCoord.from_name(config["source"]["source_name"])
    else:  # get position from ra/dec [deg]
        source_pos = SkyCoord(
            config["source"]["source_ra"],
            config["source"]["source_dec"],
            frame="icrs",
            unit="deg",
        )

    obs_ids = observations.ids

    on_region = CircleSkyRegion(center=source_pos, radius=Angle("{} deg".format(theta)))

    energy_axis = MapAxis.from_energy_bounds(
        Emin, Emax, nbin=nbins, per_decade=True, unit="TeV", name="energy"
    )

    energy_axis_true = MapAxis.from_energy_bounds(
        0.1, 20, nbin=nbins * 2, per_decade=True, unit="TeV", name="energy_true"
    )

    geom = RegionGeom.create(region=on_region, axes=[energy_axis])

    # apply exclusion regions
    exclusion_regions = []

    exclusion_regions.append(
        CircleSkyRegion(
            center=source_pos, radius=config["sky_map"]["on_exclusion_region"] * u.deg
        )
    )

    if (
        len(config["sky_map"]["exclusion_regions"]) > 0
    ):  # should be a list of CircleSkyRegions
        for region in config["sky_map"]["exclusion_regions"]:
            ra, dec = region[0]
            radius = region[1]
            exclusion_regions.append(
                CircleSkyRegion(
                    center=SkyCoord(ra, dec, unit="deg", frame="icrs"),
                    radius=radius * u.deg,
                )
            )

    # exclude nearby HAWC sources
    # hawc = SourceCatalog3HWC("$GAMMAPY_DATA/catalogs/3HWC.ecsv")
    # hawc_mask = np.sqrt(
    #        (hawc.table["ra"] - source_pos.ra.deg)**2 +
    #        (hawc.table["dec"] - source_pos.dec.deg)**2
    #        ) < 2.5
    # for src in hawc.table[hawc_mask]:
    #    if src['search_radius'] > 0:
    #        exclusion_regions.append(CircleSkyRegion(center=SkyCoord(src['ra'],src['dec'],unit='deg',frame='icrs'),radius=src['search_radius']*u.deg))
    #    else:
    #        exclusion_regions.append(CircleSkyRegion(center=SkyCoord(src['ra'],src['dec'],unit='deg',frame='icrs'),radius=0.35*u.deg))

    # exclude bright stars with 0.3 deg region (same as ED)
    star_data = np.loadtxt(
        "$GAMMAPY_DATA/catalogs/Hipparcos_MAG8_1997.dat", usecols=(0, 1, 2, 3, 4)
    )
    star_cat = Table(
        {
            "ra": star_data[:, 0],
            "dec": star_data[:, 1],
            "id": star_data[:, 2],
            "mag": star_data[:, 3],
            "colour": star_data[:, 4],
        }
    )
    star_mask = (
        np.sqrt(
            (star_cat["ra"] - source_pos.ra.deg) ** 2
            + (star_cat["dec"] - source_pos.dec.deg) ** 2
        )
        < 2.0
    )
    star_mask &= (star_cat["mag"] + star_cat["colour"]) < 8

    for src in star_cat[star_mask]:
        exclusion_regions.append(
            CircleSkyRegion(
                center=SkyCoord(src["ra"], src["dec"], unit="deg", frame="icrs"),
                radius=0.3 * u.deg,
            )
        )

    dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

    dataset_maker = SpectrumDatasetMaker(
        containment_correction=True, selection=["counts", "exposure", "edisp"]
    )

    geom = WcsGeom.create(
        npix=(150, 150), binsz=0.05, skydir=source_pos, proj="TAN", frame="icrs"
    )

    for i, exclusion_region in enumerate(exclusion_regions):
        if i == 0:
            exclusion_mask = ~geom.region_mask([exclusion_region])
        else:
            exclusion_mask *= ~geom.region_mask([exclusion_region])
    if len(exclusion_regions) > 0:
        exclusion_mask.plot()
        plt.show()
    # create reflected regions background
    # ToDo: add exclusion regions
    bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)
    safe_mask_masker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)
    datasets = Datasets()

    for obs_id, observation in zip(obs_ids, observations):
        dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
        dataset_on_off = bkg_maker.run(dataset, observation)
        dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
        datasets.append(dataset_on_off)

    info_table = datasets.info_table(cumulative=True)
    time = info_table["livetime"].to("h")
    sig = info_table["sqrt_ts"]

    # make spectrum model from user input
    if spectrum == "PL":
        # power law spectral model
        amp, idx = spectral_params
        spectral_model = PowerLawSpectralModel(
            amplitude=float(amp) * u.Unit("cm-2 s-1 TeV-1"),
            index=float(idx),
            reference=1 * u.TeV,
        )
    elif spectrum == "EXPPL":
        # exponential cutoff power law
        amp, idx, lamb = spectral_params
        spectral_model = ExpCutoffPowerLawSpectralModel(
            amplitude=float(amp) * u.Unit("cm-2 s-1 TeV-1"),
            index=float(idx),
            lambda_=float(lamb) * u.Unit("TeV-1"),
            reference=1 * u.TeV,
        )
    elif spectrum == "LP":
        amp, alp, bet = spectral_params
        spectral_model = LogParabolaSpectralModel(
            alpha=float(alp), amplitude=float(amp), reference=1 * u.TeV, beta=float(bet)
        )

    model = SkyModel(spectral_model=spectral_model, name="my_source")
    datasets.models = [model]
    fit_joint = Fit()
    result_joint = fit_joint.run(datasets=datasets)
    model_best_joint = model.copy()
    energy_edges = np.geomspace(Emin, Emax, nbins) * u.TeV

    fpe = FluxPointsEstimator(
        energy_edges=energy_edges, source="my_source", selection_optional="all"
    )
    flux_points = fpe.run(datasets=datasets)

    if plot:
        flux_points_dataset = FluxPointsDataset(
            data=flux_points, models=model_best_joint
        )
        flux_points_dataset.plot_fit()
        plt.show()
        plt.savefig(
            config["plot_names"] + "_spectrum.png", bbox_inches="tight", format="png"
        )

    return flux_points, result_joint.models, time, sig


def get_flux_lc(config, type="flux"):
    """Output run-wise flux points and the overall flux of a 1D dataset

    Parameters
    ----------
    observations: list of gammapy Observations object
    config: configuration file
    type: default = 'flux' (integral flux for whole runlist),
                    'runwise' (run by run flux points), 'custom' (time bins from config)

    Returns
    ----------
    flux_points: flux points object
    """

    theta = config["sky_map"]["theta"]
    datastore = DataStore.from_dir(config["io"]["out_dir"])
    observations = datastore.get_observations()

    amp, idx = config["spectrum"]["params"]
    if config["source"]["use_name"]:  # get position from name
        source_pos = SkyCoord.from_name(config["source"]["source_name"])
    else:  # get position from ra/dec [deg]
        source_pos = SkyCoord(
            config["source"]["source_ra"],
            config["source"]["source_dec"],
            frame="icrs",
            unit="deg",
        )
    Emin = config["spectrum"]["Emin"]
    # Emax = config["spectrum"]["Emax"]
    nbin = config["spectrum"]["nbins"]

    # selection = dict(
    #     type="sky_circle",
    #     frame="icrs",
    #     lon=source_pos.ra,
    #     lat=source_pos.dec,
    #     radius=2 * u.deg,
    # )

    # energy binning
    energy_axis = MapAxis.from_energy_bounds(
        str(Emin) + " TeV", "20 TeV", nbin
    )  # make sure all flux > Emin makes it in
    energy_axis_true = MapAxis.from_energy_bounds(
        "0.1 TeV", "30 TeV", nbin=nbin * 2, name="energy_true"
    )
    on_region_radius = Angle(str(theta) + " deg")
    on_region = CircleSkyRegion(center=source_pos, radius=on_region_radius)

    # exclusion regions
    exclusion_regions = []

    if (
        len(config["sky_map"]["exclusion_regions"]) > 0
    ):  # should be a list of CircleSkyRegions
        for region in config["sky_map"]["exclusion_regions"]:
            ra, dec = region[0]
            radius = region[1]
            exclusion_regions.append(
                CircleSkyRegion(
                    center=SkyCoord(ra, dec, unit="deg", frame="icrs"),
                    radius=radius * u.deg,
                )
            )

    # exclude nearby HAWC sources
    # hawc = SourceCatalog3HWC("$GAMMAPY_DATA/catalogs/3HWC.ecsv")
    # hawc_mask = np.sqrt(
    #        (hawc.table["ra"] - source_pos.ra.deg)**2 +
    #        (hawc.table["dec"] - source_pos.dec.deg)**2
    #        ) < 2.5
    # for src in hawc.table[hawc_mask]:
    #    if src['search_radius'] > 0:
    #        exclusion_regions.append(CircleSkyRegion(center=SkyCoord(src['ra'],src['dec'],unit='deg',frame='icrs'),radius=src['search_radius']*u.deg))
    #    else:
    #        exclusion_regions.append(CircleSkyRegion(center=SkyCoord(src['ra'],src['dec'],unit='deg',frame='icrs'),radius=0.35*u.deg))

    # exclude bright stars with 0.3 deg region (same as ED)
    star_data = np.loadtxt(
        "$GAMMAPY_DATA/catalogs/Hipparcos_MAG8_1997.dat", usecols=(0, 1, 2, 3)
    )
    star_cat = Table(
        {
            "ra": star_data[:, 0],
            "dec": star_data[:, 1],
            "id": star_data[:, 2],
            "mag": star_data[:, 3],
        }
    )
    star_mask = (
        np.sqrt(
            (star_cat["ra"] - source_pos.ra.deg) ** 2
            + (star_cat["dec"] - source_pos.dec.deg) ** 2
        )
        < 2.0
    )
    star_mask &= star_cat["mag"] < 8

    for src in star_cat[star_mask]:
        exclusion_regions.append(
            CircleSkyRegion(
                center=SkyCoord(src["ra"], src["dec"], unit="deg", frame="icrs"),
                radius=0.3 * u.deg,
            )
        )

    geom = WcsGeom.create(
        npix=(150, 150), binsz=0.05, skydir=source_pos, proj="TAN", frame="icrs"
    )

    for i, exclusion_region in enumerate(exclusion_regions):
        if i == 0:
            exclusion_mask = ~geom.region_mask([exclusion_region])
        else:
            exclusion_mask *= ~geom.region_mask([exclusion_region])

    geom = RegionGeom.create(region=on_region, axes=[energy_axis])
    dataset_maker = SpectrumDatasetMaker(
        containment_correction=True, selection=["counts", "exposure", "edisp"]
    )
    bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)
    safe_mask_masker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)

    start = Time(observations[0].gti.time_start[0], format="mjd")
    stop = Time(observations[-1].gti.time_stop[0], format="mjd")
    start.format = "isot"
    stop.format = "isot"

    if type == "flux":
        time_intervals = [Time([start, stop])]
        lc_maker_1d = LightCurveEstimator(
            energy_edges=[Emin, 20] * u.TeV,
            time_intervals=time_intervals,
            n_sigma_ul=2,
            reoptimize=False,
        )
        short_observations = observations.select_time(time_intervals)

    if type == "custom":
        binsz = config["lightcurve"]["bin_size_min"] * u.min  # bin size in minutes
        duration = (stop - start).sec / 60 * u.min
        n = int(duration / binsz)
        times = start + np.arange(n) * binsz
        time_intervals = [
            Time([tstart, tstop]) for tstart, tstop in zip(times[:-1], times[1:])
        ]
        lc_maker_1d = LightCurveEstimator(
            energy_edges=[Emin, 20] * u.TeV,
            selection_optional=None,
            time_intervals=time_intervals,
            n_sigma_ul=2,
        )
        short_observations = observations.select_time(time_intervals)

    if type == "runwise":
        lc_maker_1d = LightCurveEstimator(
            energy_edges=[Emin, 20] * u.TeV,
            selection_optional=None,
            n_sigma_ul=2,
        )
        short_observations = observations

    datasets = Datasets()

    dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

    for obs in short_observations:
        dataset = dataset_maker.run(dataset_empty.copy(), obs)
        dataset_on_off = bkg_maker.run(dataset, obs)
        dataset_on_off = safe_mask_masker.run(dataset_on_off, obs)
        datasets.append(dataset_on_off)

    spectral_model = PowerLawSpectralModel(
        index=float(idx),
        amplitude=float(amp) * u.Unit("1 / (cm2 s TeV)"),
        reference=1 * u.TeV,
    )
    spectral_model.parameters["index"].frozen = False
    spectral_model.parameters["reference"].frozen = True
    sky_model = SkyModel(
        spatial_model=None, spectral_model=spectral_model, name="model"
    )
    datasets.models = sky_model

    lc_1d = lc_maker_1d.run(datasets)
    return lc_1d
