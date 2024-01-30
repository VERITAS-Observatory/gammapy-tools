import numpy as np
from astropy.io import fits
from astropy.table import vstack
from gammapy.data import DataStore
from gammapy.irf import Background2D
from multiprocess import Pool, cpu_count

from .make_background import findData_mimic


def kl_divergence(data1, data2):
    kl = data1 * np.log(data1 / data2)
    kl = kl.sum()
    return kl


def get_background(filename):
    with fits.open(filename) as hdul:
        back = Background2D.from_hdulist(hdul)
        data = back.data

    data /= data.sum()
    data += 1e-9
    return data


def analyze_data(data, obs, sub_tab, search_dir):
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
    obs_id, config, output_name, search_runs=None, bmimic=False, overwrite=False
):
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

    with Pool(cpu_count() - 1) as pool:
        data_interest = get_background(finterest)

        def call_obs(x):
            return analyze_data(data_interest, x, sub_tab, search_dir)

        res = pool.map(call_obs, obs)

    # print (len(res))
    obs_info["kl_div"] = res
    obs_info.write(output_name, overwrite=overwrite)
    return res
