import numpy as np
from astropy.io.fits import HDUList
from astropy.time import Time
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord


def get_obs_details(hdul: HDUList) -> tuple[float, float, float]:
    """Find the run conditions

    Parameters
    ----------
    hdul (HDUList)                              - HDUL of the dl3 file

    Returns
    ----------
        zen (float)                             - Zenith angle of the observation
        nsb (float)                             - NSB Level of the observation
        az  (float)                             - Azimuth angle of the observation
    """

    zen = 90 - float(hdul[1].header["ALT_PNT"])
    nsb = float(hdul[1].header["NSBLEVEL"])
    az = float(hdul[1].header["AZ_PNT"])

    return zen, az, nsb


def get_epoch(hdul: HDUList) -> tuple[Time, float, dict]:
    """Find the epoch cooresponding to the run

    Parameters
    ----------
    hdul (HDUList)                              - HDUL of the dl3 file

    Returns
    ----------
        tobs (astropy.time.Time)                - Time of the observation
        obs_id (str)                            - Observations ID
        season (dict)                            - Relevant Epoch dictionary entry
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


def find_data_mimic(
    hdul: HDUList, config: dict, obs_table: Table
) -> tuple[np.ndarray, float]:
    """Find background data from an obs_table

    Parameters
    ----------
        hdul (HDUList)                          - HDUL of the dl3 file
        config (dict)                           - Cofigurations details
        obs_table (astropy.table.Table)         - Observations table



    Returns
    ----------
        data_mask (numpy.ndarray)               - Mask to be applied to the obs_table
        livetime (float)                        - Livetime of the selected data

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
    sub_table = obs_table[np.abs(obs_date - tobs) < time_max]

    # mask in NSB
    sub_table = sub_table[np.abs(sub_table["NSBLEVEL"] - nsb_obs) < nsb_diff]
    # mask in zenith
    el = np.deg2rad(90 - sub_table["ZEN_PNT"])
    sub_table = sub_table[np.abs(1 / np.sin(el) - 1 / np.sin(el_obs)) < el_diff]
    # Mask in az
    sub_table = sub_table[np.abs(sub_table["AZ_PNT"] - az_obs) < az_diff]
    # Mask out < 4 tels
    sub_table = sub_table[sub_table["N_TELS"] == n_tel]
    # Mask out galactic plane
    sub_table = sub_table[
        (
            np.abs(SkyCoord(sub_table["RA_PNT"], sub_table["DEC_PNT"]).galactic.b)
            > 10 * u.deg
        )
    ]
    # Make sure we're not using the same run!
    sub_table = sub_table[sub_table["OBS_ID"] != obs_id]

    # Cutting on elevation
    if "el_min" in config["background_selection"]:
        el = np.deg2rad(90 - sub_table["ZEN_PNT"])

        sub_table = sub_table[el > config["background_selection"]["el_min"]]
    if "el_max" in config["background_selection"]:
        el = np.deg2rad(90 - sub_table["ZEN_PNT"])

        sub_table = sub_table[el < config["background_selection"]["el_max"]]

    livetime = np.sum(sub_table["LIVETIME"] / 60 / 60)
    obs_list = sub_table["OBS_ID"]
    return obs_list, livetime
