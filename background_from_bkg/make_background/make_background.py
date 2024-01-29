import os
import sys
from multiprocessing import Pool
from os import listdir
from os.path import isfile, join

import numpy as np
import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord

# Astropy
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.time import Time
from gammapy.data import DataStore
from gammapy.irf import Background2D

# Gammapy stuff
from gammapy.maps import MapAxis

# V2DL3 stuff
from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file

from .background_models import BackgroundModelEstimator


def get_obs_details(hdul):
    """Find the run conditions

    Parameters
    ----------
    hdul : HDUList of the dl3 file

    Returns
    ----------
        zen  - Zenith angle of the observation
        nsb   - NSB Level of the observation
        az    - Azimuth angle of the observation
    """

    zen = 90 - float(hdul[1].header["ALT_PNT"])
    nsb = float(hdul[1].header["NSBLEVEL"])
    az = float(hdul[1].header["AZ_PNT"])

    return zen, az, nsb


def get_epoch(hdul) -> (float, list):
    """Find the epoch cooresponding to the run

    Parameters
    ----------
    hdul : HDUList of the dl3 file

    Returns
    ----------
        tobs    - Time of the observation
        obs_id  - Observations ID
        season  - Relevant Epoch dictionary entry
    """
    epochs = {
        "V4": {"tstart": Time("2000-01-01"), "tstop": Time("2009-09-13")},
        "V5": {"tstart": Time("2009-09-14"), "tstop": Time("2012-07-31")},
        "V6_2012_2013a": {"tstart": Time("2012-08-01"), "tstop": Time("2013-03-15")},
        "V6_2012_2013b": {"tstart": Time("2013-03-16"), "tstop": Time("2013-11-16")},
        "V6_2013_2014a": {"tstart": Time("2013-11-17"), "tstop": Time("2014-05-12")},
        "V6_2013_2014b": {"tstart": Time("2014-05-13"), "tstop": Time("2014-11-07")},
        "V6_2014_2015": {"tstart": Time("2014-11-08"), "tstop": Time("2015-07-31")},
        "V6_2015_2016": {"tstart": Time("2015-08-01"), "tstop": Time("2016-07-31")},
        "V6_2016_2017": {"tstart": Time("2016-08-01"), "tstop": Time("2017-07-31")},
        "V6_2017_2018": {"tstart": Time("2017-08-01"), "tstop": Time("2018-07-31")},
        "V6_2018_2019": {"tstart": Time("2018-08-01"), "tstop": Time("2019-11-12")},
        "V6_2019_2020w": {"tstart": Time("2019-11-13"), "tstop": Time("2020-05-07")},
        "V6_2020_2020s": {"tstart": Time("2020-05-08"), "tstop": Time("2020-11-04")},
        "V6_2020_2021w": {"tstart": Time("2020-11-05"), "tstop": Time("2021-04-27")},
        "V6_2021_2021s": {"tstart": Time("2021-04-28"), "tstop": Time("2021-11-16")},
        "V6_2021_2022w": {"tstart": Time("2021-11-17"), "tstop": Time("2022-05-09")},
        #  "V6_2022_2022s" : {"tstart"  :Time("2022-05-10"), "tstop" : Time("2022-11-08")},
        "V6_2022_2022s": {"tstart": Time("2022-05-10"), "tstop": Time("2025-11-08")},
    }

    tobs = Time(hdul[1].header["DATE-OBS"])
    # for e ine.keys():
    #     if e[e]["tstart"] < tobs) and e[e]["tstop"] > tobs):
    #         returne[e]["tstart"],e[e]["tstop"], tobs
    obs_id = hdul[1].header["OBS_ID"]
    season = ""

    for key, value in epochs.items():
        if (tobs.iso > value["tstart"]) & (tobs.iso <= value["tstop"]):
            season = epochs[key]
    # to do: add error catching in case the time specified is outside of the VERITAS epochs
    return tobs, obs_id, season


def find_data_mimic(hdul, config, obs_table) -> (list, float):
    """Find background data from an obs_table

    Parameters
    ----------
        hdul    - HDUList of the dl3 file
        config  - Dictionary containing cofigurations details



    Returns
    ----------
        data_mask   - Mask to be applied to the obs_table with selected data
        livetime    - Livetime of the selected data

    """
    try:
        tobs, obs_id, season = get_epoch(hdul)

    # Should only happen if there is an issue with the dl3 file...
    # Todo: better handle here
    except Exception as e:
        print(e)
        print(hdul)
        print(hdul[0].header)
        return None

    # obs_table = DataStore.from_dir(config["io"]["search_datastore"]).obs_table

    zen_obs, az_obs, nsb_obs = get_obs_details(hdul)
    el_obs = np.deg2rad(90 - zen_obs)

    obs_date = Time(obs_table["DATE-AVG"])
    el = np.deg2rad(90 - obs_table["ZEN_PNT"])

    # Getting parameters or setting the defaults
    el_diff = (
        1.0 / np.sin(config["background_selection"]["el_diff"])
        if "el_diff" in config["background_selection"]
        else 0.15
    )
    az_diff = (
        config["background_selection"]["az_diff"]
        if "az_diff" in config["background_selection"]
        else 30.0
    )
    nsb_diff = (
        config["background_selection"]["nsb_diff"]
        if "nsb_diff" in config["background_selection"]
        else 0.7
    )
    time_max = (
        config["background_selection"]["time_max"]
        if "time_max" in config["background_selection"]
        else 30
    )
    n_tel = (
        config["background_selection"]["n_tel"]
        if "n_tel" in config["background_selection"]
        else 4
    )

    # mask in time
    data_mask = np.abs(obs_date - tobs) < time_max

    # mask in NSB
    data_mask &= np.abs(obs_table["NSBLEVEL"] - nsb_obs) < nsb_diff
    # mask in zenith
    data_mask &= np.abs(1 / np.sin(el) - 1 / np.sin(el_obs)) < el_diff
    # Mask in az
    data_mask &= np.abs(obs_table["AZ_PNT"] - az_obs) < az_diff
    # Mask out < 4 tels
    data_mask &= obs_table["N_TELS"] == n_tel
    # Mask out galactic plane
    data_mask &= (
        np.abs(SkyCoord(obs_table["RA_PNT"], obs_table["DEC_PNT"]).galactic.b)
        > 10 * u.deg
    )
    # Make sure we're not using the same run!
    data_mask &= obs_table["OBS_ID"] != obs_id

    # Cutting on elevation
    if "el_min" in config["background_selection"]:
        data_mask &= el > config["background_selection"]["el_min"]
    if "el_max" in config["background_selection"]:
        data_mask &= el < config["background_selection"]["el_max"]

    livetime = np.sum(obs_table[data_mask]["LIVETIME"] / 60 / 60)

    return data_mask, livetime


# ToDo: are table operations cpu bound? mega_store = path / mega_store = table
# Are tables "pickleable?"
def get_background_for_run(parms) -> (str, list):
    """Generate the background for a given run

    Parameters
    ----------
    parms   - tupple of obs_id and config
        obs_id - Run number
        config  - Dictionary containing cofigurations details


    Returns
    ----------
        out_file    - Filename of the output file
        obs_list    - List of obseravtions used to generate the background

    """

    obs, config = parms

    in_dir = config["io"]["in_dir"]
    out_dir = config["io"]["out_dir"]

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    mega_store = DataStore.from_dir(config["io"]["search_datastore"])
    mega_table = mega_store.obs_table
    local_store = DataStore.from_dir(in_dir)

    tab_local = local_store.hdu_table
    mask = tab_local["OBS_ID"] == obs
    file_name = str(tab_local[mask]["FILE_NAME"][0])
    file_path = str(tab_local[mask]["FILE_DIR"][0])

    in_file = os.path.join(in_dir, file_path, file_name)
    # print (in_file, file_path, file_name)
    with fits.open(in_file) as hdul:
        # try/catch should only trigger if run is outside of epoch definitions
        # or we have a major issue with the DL3 file
        # ToDo: Better error exiting
        # try :

        if "bkg_runlist" not in config["background_selection"]:
            config["background_selection"]["bkg_runlist"] = {}

        if obs not in config["background_selection"]["bkg_runlist"]:
            data_mask, livetime = find_data_mimic(hdul, config, mega_table)
            if livetime < 10:
                print(obs, livetime)
            obs_list = mega_table[data_mask]["OBS_ID"]
            config["background_selection"]["bkg_runlist"][obs] = obs_list

        else:
            obs_list = config["background_selection"]["bkg_runlist"][obs]

        # Get observations form the mega store
        observations = mega_store.get_observations(
            obs_list, required_irf="all-optional"
        )

        # Defined energy axes
        energy = MapAxis.from_energy_bounds(
            config["binning"]["e_min"],
            config["binning"]["e_max"],
            config["binning"]["e_bins"],
            name="energy",
            unit="TeV",
        )

        offset = MapAxis.from_bounds(
            config["binning"]["off_min"],
            config["binning"]["off_max"],
            nbin=config["binning"]["off_bins"],
            interp="lin",
            unit="deg",
            name="offset",
        )

        estimator = BackgroundModelEstimator(
            energy, offset, smooth=config["background_selection"]["smooth"]
        )
        estimator.run(observations)

        hdul.append(estimator.background_rate.to_table_hdu())

        out_file = os.path.join(out_dir, file_name)

        hdul.writeto(out_file, overwrite=True)

    return out_file, obs_list


def generate_background_from_run(parms) -> str:
    """Generate the background from a given run

    Parameters
    ----------
    parms   - tupple of obs_id and config
        obs_id - Run number
        config  - Dictionary containing cofigurations details


    Returns
    ----------
        data_mask   - Mask to be applied to the obs_table with selected data
        livetime    - Livetime of the selected data

    """
    obs, config = parms

    in_dir = config["io"]["in_dir"]
    out_dir = config["io"]["out_dir"]

    # Assume a datastore
    data_store = DataStore.from_dir(config["io"]["in_dir"])
    # data_table = data_store.obs_table

    tab = data_store.hdu_table
    mask = tab["OBS_ID"] == obs
    file_name = str(tab[mask]["FILE_NAME"][0])
    file_path = str(tab[mask]["FILE_DIR"][0])

    in_file = os.path.join(in_dir, file_path, file_name)

    with fits.open(in_file) as hdul:
        # try/catch should only trigger if run is outside of epoch definitions
        # or we have a major issue withe the DL3 file
        # ToDo: Better error exiting
        # try :
        # except Exception as e:
        #     print (e)
        #     return None

        # Get observations form the mega store
        observations = data_store.get_observations([obs], required_irf="all-optional")

        # Defined energy axes
        energy = MapAxis.from_energy_bounds(
            config["binning"]["e_min"],
            config["binning"]["e_max"],
            config["binning"]["e_bins"],
            name="energy",
            unit="TeV",
        )

        offset = MapAxis.from_bounds(
            config["binning"]["off_min"],
            config["binning"]["off_max"],
            nbin=config["binning"]["off_bins"],
            interp="lin",
            unit="deg",
            name="offset",
        )

        estimator = BackgroundModelEstimator(
            energy, offset, smooth=config["background_selection"]["smooth"]
        )
        estimator.run(observations)

        hdul.append(estimator.background_rate.to_table_hdu())
        out_file = os.path.join(out_dir, file_name)

        hdul.writeto(out_file, overwrite=True)
    return out_file


def write_index_files(config):
    dl3_dir = config["io"]["out_dir"]
    dl3Files = [
        dl3_dir + f
        for f in listdir(dl3_dir)
        if isfile(join(dl3_dir, f))
        and f.endswith(".fits")
        and (f.strip(".anasum.fits"))
    ]
    create_obs_hdu_index_file(dl3Files, index_file_dir=dl3_dir)
    return


def attach_bkg(parms):
    # this function is only used to append background files you already have
    # (i.e., not from the mega store) - not sure how useful this is??
    obs_id, config = parms
    bkg_dir = config["io"]["out_dir"]
    dl3_dir = config["io"]["in_dir"]

    dl3Files = [
        dl3_dir + f
        for f in listdir(dl3_dir)
        if isfile(join(dl3_dir, f))
        and f.endswith(".fits")
        and (f.strip(".anasum.fits"))
    ]
    create_obs_hdu_index_file(dl3Files, index_file_dir=dl3_dir)
    data_store = DataStore.from_dir(dl3_dir)
    hdu_table = data_store.hdu_table.read(f"{dl3_dir}" + "/" + "hdu-index.fits.gz")
    obs_table = data_store.obs_table.read(f"{dl3_dir}" + "/" + "obs-index.fits.gz")
    data_store.hdu_table.remove_rows(data_store.hdu_table["HDU_TYPE"] == "bkg")
    hdu_table.remove_rows(hdu_table["HDU_NAME"] == "BKG")
    bkg_model_names = "background_" + str(obs_id) + "_smoothed.fits"

    if "bkg" not in list(hdu_table["HDU_TYPE"]):
        hdu_table.remove_column("FILE_DIR")
        hdu_table.add_column(f"{dl3_dir}", name="FILE_DIR")
        # zenith_angle = list(obs_table["ZEN_PNT"])
        # azimuth_angle = list(obs_table["AZ_PNT"])
        # obs_ids = list(obs_table["OBS_ID"])
        rows = []
        filename = bkg_model_names
        row = []
        row = {
            "OBS_ID": obs_id,
            "HDU_TYPE": "bkg",
            "HDU_CLASS": "bkg_2d",
            "FILE_DIR": f"{bkg_dir}",
            "FILE_NAME": filename,
            "HDU_NAME": "BKG",
        }
        rows.append(row)
        hdu_table_bkg = Table(rows=rows)
        hdu_table = vstack([hdu_table, hdu_table_bkg])
        hdu_table.sort("OBS_ID")
        hdu_bin_table = fits.table_to_hdu(hdu_table)
        hdul = fits.open(f"{dl3_dir}" + "/" + "hdu-index.fits.gz")
        hdul[1] = hdu_bin_table
        hdul.writeto(f"{dl3_dir}" + "/" + "hdu-index.fits.gz", overwrite=True)
    filename = list(hdu_table["FILE_NAME"])
    file_dir = list(hdu_table["FILE_DIR"])
    bkg_indl = []
    for i in range(len(file_dir)):
        if "background" in filename[i]:
            bkg_indl.append(i)

    delete_prior_bkg = True
    count = 0

    hdu_list = fits.open(f"{dl3_dir}" + "/" + str(obs_id) + ".anasum.fits")

    hdul_new = fits.HDUList()
    for index, hdu in enumerate(hdu_list):
        # remove BKG and replace by new one
        if index > 0:  # first index not containing EXTNAME
            if delete_prior_bkg and hdu.header["EXTNAME"] == "BKG":
                print(index, "removed")
                continue
        # hdu_copy = hdu.copy()
        hdul_new.append(hdu)
    bkg_index = bkg_indl[count]
    print(hdu_table["FILE_DIR"][bkg_index] + hdu_table["FILE_NAME"][bkg_index])
    BKG2D = Background2D.read(
        hdu_table["FILE_DIR"][bkg_index] + hdu_table["FILE_NAME"][bkg_index]
    )
    hdu = BKG2D.to_table_hdu()
    hdu.name = "BKG"
    hdul_new.append(hdu)
    count += 1
    print([_.name for _ in hdul_new])
    hdul_new.writeto(f"{bkg_dir}" + "/" + str(obs_id) + ".bkg.fits", overwrite=True)
    data, header = fits.getdata(
        f"{bkg_dir}" + "/" + str(obs_id) + ".bkg.fits", 6, header=True
    )
    header["HDUCLAS1"] = "RESPONSE"
    header["HDUCLAS2"] = "BKG"
    header["HDUCLAS3"] = "FULL-ENCLOSURE"
    header["HDUCLAS4"] = "BKG_2D"
    fits.update(f"{bkg_dir}" + "/" + str(obs_id) + ".bkg.fits", data, header, 6)
    final = []

    for i in list(obs_table["OBS_ID"]):
        a = f"{bkg_dir}" "/" + str(i) + ".bkg.fits"
        final.append(a)

    create_obs_hdu_index_file(final, index_file_dir=bkg_dir)
    return ()


def get_mimic_for_run(parms):
    obs, config = parms

    dl3_file_fmt = "{runid}.anasum.fits"

    in_dir = config["io"]["in_dir"]
    # out_dir = config["io"]["out_dir"]

    mega_store = DataStore.from_dir(config["io"]["search_datastore"])
    mega_table = mega_store.obs_table

    in_file = os.path.join(in_dir, dl3_file_fmt.format(runid=obs))
    with fits.open(in_file) as hdul:
        data_mask, livetime = find_data_mimic(hdul, config)

        if livetime < 10:
            print(obs, livetime)

        # Get observations form the mega store
        observations = mega_store.get_observations(
            mega_table[data_mask]["OBS_ID"], required_irf="all-optional"
        )
    return observations


def run_make_background(config) -> dict:
    """Generate background

    Parameters
    ----------
    config   - dict
        dictionary containing cofigurations details



    Returns
    ----------
    config   - dict
        dictionary containing cofigurations details

    """
    in_dir = config["io"]["in_dir"]
    out_dir = config["io"]["out_dir"]
    from_run = config["io"]["from_run"]

    # Create a directory if it doesn't exist
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Targeted data store
    data_store = DataStore.from_dir(in_dir)
    obs_table = data_store.obs_table

    flist = []

    # if "bkg_runlist" in config["background_selection"]:
    #     obs_parameters =  [(obs, config) for obs in config["background_selection"]["bkg_runlist"]]
    # else:
    obs_parameters = [(obs, config) for obs in obs_table["OBS_ID"]]

    if from_run:
        map_fn = generate_background_from_run
    else:
        map_fn = get_background_for_run

    with Pool(config["config"]["njobs"]) as pool:
        output = pool.map(map_fn, obs_parameters)

    if from_run:
        flist = output
    else:
        flist = [fl[0] for fl in output if fl is not None]

    create_obs_hdu_index_file(flist, out_dir)
    config["background_selection"]["bkg_runlist"] = {}
    for obs, fl in zip(list(obs_table["OBS_ID"]), output):
        # print (obs, fl[0])
        # print (fl[1])
        config["background_selection"]["bkg_runlist"][obs] = list(fl[1])

    return config


if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        config = yaml.safe_load(f)

    run_make_background(config)
