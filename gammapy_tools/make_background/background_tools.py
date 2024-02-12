import numpy as np
from astropy.io import fits
from astropy.table import vstack, Table
from gammapy.data import DataStore
from gammapy.irf import Background2D
from multiprocess import Pool, cpu_count

from .make_background import findData_mimic


def kl_divergence(data1: np.array, data2: np.array) -> float:
    """Calculate the Kullback-Leibler Divergence.
    This provides a metric for comparing two 2D distributions

    Parameters
    ----------
        data1 (numpy.array)         - Array to compare against.
        data2 (numpy.array)         - Array to compare.


    Returns
    ----------
        kl  (float)                 - KL Divergence score.

    """
    kl = data1 * np.log(data1 / data2)
    kl = kl.sum()
    return kl


def get_background(filename: str) -> np.array:
    """Get the background from a fits file

    Parameters
    ----------
        filename (str)              - Name of the file to be read in.


    Returns
    ----------
        data (numpy.array)          - Background.

    """
    with fits.open(filename) as hdul:
        back = Background2D.from_hdulist(hdul)
        data = back.data

    data /= data.sum()
    data += 1e-9
    return data


def analyze_data(data: np.array, obs: int, sub_tab: Table, search_dir: str) -> float:
    """Get the KL Divergence for an observation

    Parameters
    ----------
        data (numpy.array)          - Reference distribution.
        obs (int)                   - Observation ID of interest.
        sub_tab (astropy.table)     - Table of observations for searching.
        search_dir (str)            - Directory to search for the observation.


    Returns
    ----------
        kl_div (float)              - KL Divergence.

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
    output_name: str,
    search_runs: list = None,
    bmimic: bool = False,
    overwrite: bool = False,
    ncpu: int = None,
) -> list:
    """Calculate the KL Divergence for a set of observations

    Parameters
    ----------
        obs_id (int)                - Observation of interest.
        config (dict)               - Configuration dictionary.
        output_name (str)           - Name of the file to save the results to.
        search_runs (list)          - List of runs to calculate the KL divergence for.
                                      Defaults to None
        bmimic (bool)               - If the mimic critera is used to find the `search_runs`.
                                      Defaults to False.
        overwrite (bool)            - Overwrite the file specified by `output_name`.
                                      Defaults to False.
        ncpu (int)                  - Number of parallel jobs to run. Defaults to `ncpu-1`.


    Returns
    ----------
        res (list)                  - List of KL divergence

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
            data_mask, _ = findData_mimic(hdul, config, obs_info)
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

    if ncpu is None:
        ncpu = cpu_count() - 1

    with Pool(ncpu) as pool:
        data_interest = get_background(finterest)

        def call_obs(x):
            return analyze_data(data_interest, x, sub_tab, search_dir)

        res = pool.map(call_obs, obs)

    obs_info["kl_div"] = res
    obs_info.write(output_name, overwrite=overwrite)
    return res
