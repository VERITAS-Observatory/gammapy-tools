"""
Background tools

A set of tools to find and compare backgrounds to determine the most suitable runs
for background generation
"""

import numpy as np
from astropy.io import fits
from astropy.table import vstack, Table
from gammapy.data import DataStore
from gammapy.irf import Background2D
from multiprocess import Pool, cpu_count
from ..utils import get_cdf

from ..utils.run_details import find_data_mimic


def kl_divergence(data1: np.ndarray, data2: np.ndarray) -> float:
    """Calculate the Kullback-Leibler Divergence.
    This provides a metric for comparing two 2D distributions

    Parameters
    ----------
        data1 (numpy.ndarray)                   - Array to compare against.
        data2 (numpy.ndarray)                   - Array to compare.


    Returns
    ----------
        kl  (float)                             - KL Divergence score.

    """
    kl = data1 * np.log(data1 / data2)
    kl = kl.sum()
    return kl


def get_background(filename: str) -> np.ndarray:
    """Get the background from a fits file

    Parameters
    ----------
        filename (str)                          - Name of the file to be read in.


    Returns
    ----------
        data (numpy.ndarray)                    - Background.

    """
    with fits.open(filename) as hdul:
        back = Background2D.from_hdulist(hdul)
        data = back.data

    data /= data.sum()
    data += 1e-9
    return data


def analyze_data(data: np.ndarray, obs: int, sub_tab: Table, search_dir: str) -> float:
    """Get the KL Divergence for an observation

    Parameters
    ----------
        data (numpy.ndarray)                    - Reference distribution.
        obs (int)                               - Observation ID of interest.
        sub_tab (astropy.table)                 - Table of observations for searching.
        search_dir (str)                        - Directory to search for the observation.


    Returns
    ----------
        kl_div (float)                          - KL Divergence.

    """
    of_interest = sub_tab["OBS_ID"] == obs
    fname = (
        search_dir
        + "/"
        + sub_tab[of_interest]["FILE_DIR"][0]
        + "/"
        + sub_tab[of_interest]["FILE_NAME"][0]
    )
    tmp_back = get_background(fname)
    kl_div = kl_divergence(data, tmp_back)
    return kl_div


def process_run(
    obs_id: int,
    config: dict,
    output_name: str = None,
    search_runs: list = None,
    bmimic: bool = False,
    overwrite: bool = False,
    njobs: int = None,
) -> Table:
    """Calculate the KL Divergence for a set of observations

    Parameters
    ----------
        obs_id (int)                            - Observation of interest.
        config (dict)                           - Configuration dictionary.
        output_name (str)                       - Name of the file to save the results to.
                                                  Defaults to None, if None not file is written
        search_runs (list)                      - List of runs to calculate the KL divergence for.
                                                  Defaults to None
        bmimic (bool)                           - If the mimic critera is used to find
                                                  the `search_runs`.
                                                  Defaults to False.
        overwrite (bool)                        - Overwrite the file specified by `output_name`.
                                                  Defaults to False.
        njobs (int)                             - Number of parallel jobs to run.
                                                  Defaults to `njobs-1`.


    Returns
    ----------
        obs_info (astropy.table.Table)          - Table including KL divergence

    """
    search_dir = config["io"]["search_datastore"]

    data_store = DataStore.from_dir(search_dir)

    sub_tab = data_store.hdu_table[data_store.hdu_table["HDU_CLASS"] == "bkg_2d"]
    obs_info = data_store.obs_table

    of_interest = data_store.hdu_table["OBS_ID"] == obs_id

    finterest = (
        search_dir
        + "/"
        + data_store.hdu_table[of_interest]["FILE_DIR"][0]
        + "/"
        + data_store.hdu_table[of_interest]["FILE_NAME"][0]
    )

    if search_runs is None:
        obs = data_store.obs_ids
    elif bmimic:
        with fits.open(finterest) as hdul:
            data_mask, _ = find_data_mimic(hdul, config, obs_info)
        obs = obs_info[data_mask]["OBS_ID"]
    else:
        obs = search_runs
    stack = []
    info_tab = []
    for o in obs:
        info_tab.append(obs_info[obs_info["OBS_ID"] == o])
        stack.append(sub_tab[sub_tab["OBS_ID"] == o])
    sub_tab = vstack(stack)
    obs_info = vstack(info_tab)

    if njobs is None:
        njobs = cpu_count() - 1

    with Pool(njobs) as pool:
        data_interest = get_background(finterest)

        def call_obs(x):
            return analyze_data(data_interest, x, sub_tab, search_dir)

        res = pool.map(call_obs, obs)

    obs_info["KL_DIV"] = res
    if output_name is not None:
        obs_info.write(output_name, overwrite=overwrite)
    return obs_info


def get_requested_exposure(obs_table: Table, tobs: float) -> Table:
    """Get a list of observations which match the required observations time

    Parameters
    ----------
        obs_table (astropy.table.Table)         - Table of observations.
        tobs (float)                            - Requested background exposure (hours)


    Returns
    ----------
        obs (astropy.table.Table)               - Table of runs to use

    """

    # Sort by kl_div (reverse, kl should be small)
    obs_table.sort(["KL_DIV"], reverse=False)
    # The 0th item should be the same run
    obs_table = obs_table[1:]
    # Get the cumulative livetime
    cumul_obs = get_cdf(obs_table["LIVETIME"], normalize=False)
    # Convert to hours and 1% buffer
    mask = cumul_obs / 60 / 60 <= 1.01 * tobs

    return obs_table[mask]
