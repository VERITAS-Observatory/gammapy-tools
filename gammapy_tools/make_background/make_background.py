import os
import sys
from multiprocessing import Pool

import numpy as np
import yaml
from os import listdir
from os.path import isfile, join, getsize
# Astropy
from astropy.io import fits
from gammapy.data import DataStore

# Gammapy stuff
from gammapy.maps import MapAxis

# V2DL3 stuff
from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file

# From here
from .background_models import BackgroundModelEstimator
from .background_tools import process_run, get_requested_exposure
from ..utils.run_details import find_data_mimic


# ToDo: are table operations cpu bound? mega_store = path / mega_store = table
# Are tables "pickleable?"
def get_background_for_run(parms: tuple[float, dict]) -> tuple[str, list]:
    """Generate the background for a given run

    Parameters
    ----------
        parms (tuple[float,dict])               - tupple of Run number and config



    Returns
    ----------
        out_file (str)                          - Name of the output file
        obs_list (list)                         - Obseravtions used to generate the background

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

    with fits.open(in_file) as hdul:

        if "bkg_runlist" not in config["background_selection"]:
            config["background_selection"]["bkg_runlist"] = {}

        if obs not in config["background_selection"]["bkg_runlist"]:

            if config["background_selection"]["KL_DIV"]:
                kl_table = process_run(
                    obs,
                    config,
                    output_name=None,
                    search_runs=None,
                    bmimic=True,
                    overwrite=False,
                    njobs=config["config"]["njobs"],
                )

                # Default to 10 hours
                t_req = 10
                if "time_request" in config["background_selection"]:
                    t_req = float(config["background_selection"]["time_request"])

                kl_table = get_requested_exposure(kl_table, t_req)
                livetime = np.sum(kl_table["LIVETIME"] / 60 / 60)

                if livetime < t_req:
                    print(f"{obs} used {livetime} hours for background generation")

                obs_list = kl_table["OBS_ID"]

            else:
                obs_list, livetime = find_data_mimic(hdul, config, mega_table)
                # obs_list = mega_table[data_mask]["OBS_ID"]

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

        smooth_sigma = 1
        if "smooth_sigma" in config["background_selection"]:
            smooth_sigma = config["background_selection"]["smooth_sigma"]

        estimator = BackgroundModelEstimator(
            energy,
            offset,
            smooth=config["background_selection"]["smooth"],
            smooth_sigma=smooth_sigma,
            njobs=config["config"]["njobs"],
        )

        estimator.run(observations)
        # Check if a background currently exists
        if "BACKGROUND" in hdul:
            bkg_indx = hdul.index_of("BACKGROUND")
            hdul.pop(bkg_indx)

        hdul.append(estimator.background_rate.to_table_hdu())

        out_file = os.path.join(out_dir, file_name)

        hdul.writeto(out_file, overwrite=True)

    return out_file, obs_list


def generate_background_from_run(parms: tuple[int, dict]) -> str:
    """Generate the background from a given run

    Parameters
    ----------
        parms (tuple[int, dict])                    - tupple of run number and config dictionary


    Returns
    ----------
        filename (str)                              - Name of the background attached file

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
        if "BACKGROUND" in hdul:
            bkg_indx = hdul.index_of("BACKGROUND")
            hdul.pop(bkg_indx)

        hdul.append(estimator.background_rate.to_table_hdu())
        out_file = os.path.join(out_dir, file_name)

        hdul.writeto(out_file, overwrite=True)

    return out_file


# SOB Repeated code?

# def write_index_files(config):
#     dl3_dir = config["io"]["out_dir"]
#     dl3Files = [
#         dl3_dir + f
#         for f in listdir(dl3_dir)
#         if isfile(join(dl3_dir, f))
#         and f.endswith(".fits")
#         and (f.strip(".anasum.fits"))
#     ]
#     create_obs_hdu_index_file(dl3Files, index_file_dir=dl3_dir)
#     return


# def attach_bkg(parms):
#     # this function is only used to append background files you already have
#     # (i.e., not from the mega store) - not sure how useful this is??
#     obs_id, config = parms
#     bkg_dir = config["io"]["out_dir"]
#     dl3_dir = config["io"]["in_dir"]

#     dl3Files = [
#         dl3_dir + f
#         for f in listdir(dl3_dir)
#         if isfile(join(dl3_dir, f))
#         and f.endswith(".fits")
#         and (f.strip(".anasum.fits"))
#     ]
#     create_obs_hdu_index_file(dl3Files, index_file_dir=dl3_dir)
#     data_store = DataStore.from_dir(dl3_dir)
#     hdu_table = data_store.hdu_table.read(f"{dl3_dir}" + "/" + "hdu-index.fits.gz")
#     obs_table = data_store.obs_table.read(f"{dl3_dir}" + "/" + "obs-index.fits.gz")
#     data_store.hdu_table.remove_rows(data_store.hdu_table["HDU_TYPE"] == "bkg")
#     hdu_table.remove_rows(hdu_table["HDU_NAME"] == "BKG")
#     bkg_model_names = "background_" + str(obs_id) + "_smoothed.fits"

#     if "bkg" not in list(hdu_table["HDU_TYPE"]):
#         hdu_table.remove_column("FILE_DIR")
#         hdu_table.add_column(f"{dl3_dir}", name="FILE_DIR")
#         # zenith_angle = list(obs_table["ZEN_PNT"])
#         # azimuth_angle = list(obs_table["AZ_PNT"])
#         # obs_ids = list(obs_table["OBS_ID"])
#         rows = []
#         filename = bkg_model_names
#         row = []
#         row = {
#             "OBS_ID": obs_id,
#             "HDU_TYPE": "bkg",
#             "HDU_CLASS": "bkg_2d",
#             "FILE_DIR": f"{bkg_dir}",
#             "FILE_NAME": filename,
#             "HDU_NAME": "BKG",
#         }
#         rows.append(row)
#         hdu_table_bkg = Table(rows=rows)
#         hdu_table = vstack([hdu_table, hdu_table_bkg])
#         hdu_table.sort("OBS_ID")
#         hdu_bin_table = fits.table_to_hdu(hdu_table)
#         hdul = fits.open(f"{dl3_dir}" + "/" + "hdu-index.fits.gz")
#         hdul[1] = hdu_bin_table
#         hdul.writeto(f"{dl3_dir}" + "/" + "hdu-index.fits.gz", overwrite=True)
#     filename = list(hdu_table["FILE_NAME"])
#     file_dir = list(hdu_table["FILE_DIR"])
#     bkg_indl = []
#     for i in range(len(file_dir)):
#         if "background" in filename[i]:
#             bkg_indl.append(i)

#     delete_prior_bkg = True
#     count = 0

#     hdu_list = fits.open(f"{dl3_dir}" + "/" + str(obs_id) + ".anasum.fits")

#     hdul_new = fits.HDUList()
#     for index, hdu in enumerate(hdu_list):
#         # remove BKG and replace by new one
#         if index > 0:  # first index not containing EXTNAME
#             if delete_prior_bkg and hdu.header["EXTNAME"] == "BKG":
#                 print(index, "removed")
#                 continue
#         # hdu_copy = hdu.copy()
#         hdul_new.append(hdu)
#     bkg_index = bkg_indl[count]
#     print(hdu_table["FILE_DIR"][bkg_index] + hdu_table["FILE_NAME"][bkg_index])
#     BKG2D = Background2D.read(
#         hdu_table["FILE_DIR"][bkg_index] + hdu_table["FILE_NAME"][bkg_index]
#     )
#     hdu = BKG2D.to_table_hdu()
#     hdu.name = "BKG"
#     hdul_new.append(hdu)
#     count += 1
#     print([_.name for _ in hdul_new])
#     hdul_new.writeto(f"{bkg_dir}" + "/" + str(obs_id) + ".bkg.fits", overwrite=True)
#     data, header = fits.getdata(
#         f"{bkg_dir}" + "/" + str(obs_id) + ".bkg.fits", 6, header=True
#     )
#     header["HDUCLAS1"] = "RESPONSE"
#     header["HDUCLAS2"] = "BKG"
#     header["HDUCLAS3"] = "FULL-ENCLOSURE"
#     header["HDUCLAS4"] = "BKG_2D"
#     fits.update(f"{bkg_dir}" + "/" + str(obs_id) + ".bkg.fits", data, header, 6)
#     final = []

#     for i in list(obs_table["OBS_ID"]):
#         a = f"{bkg_dir}" "/" + str(i) + ".bkg.fits"
#         final.append(a)

#     create_obs_hdu_index_file(final, index_file_dir=bkg_dir)
#     return ()


# def get_mimic_for_run(parms):
#     obs, config = parms

#     dl3_file_fmt = "{runid}.anasum.fits"

#     in_dir = config["io"]["in_dir"]
#     # out_dir = config["io"]["out_dir"]

#     mega_store = DataStore.from_dir(config["io"]["search_datastore"])
#     mega_table = mega_store.obs_table

#     in_file = os.path.join(in_dir, dl3_file_fmt.format(runid=obs))
#     with fits.open(in_file) as hdul:
#         data_mask, livetime = find_data_mimic(hdul, config)

#         if livetime < 10:
#             print(obs, livetime)

#         # Get observations form the mega store
#         observations = mega_store.get_observations(
#             mega_table[data_mask]["OBS_ID"], required_irf="all-optional"
#         )
#     return observations


def run_make_background(config: dict) -> dict:
    """Generate background

    Parameters
    ----------
        config (dict)                       - Configuration details

    Returns
    ----------
        config (dict)                       - Modified configuration of background attached runs
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
        with Pool(config["config"]["njobs"]) as pool:
            output = pool.map(map_fn, obs_parameters)
    else:
        map_fn = get_background_for_run
        with Pool(1) as pool:
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

def write_index_files(config):
    dl3_dir = config["io"]["out_dir"]
    dl3Files = [dl3_dir + f for f in listdir(dl3_dir) if isfile(join(dl3_dir, f)) and (f.endswith(".fits") or (f.endswith(".fits.gz") and not f.startswith("obs") and not f.startswith("hdu"))) and (f.strip('.anasum.fits') )]
    create_obs_hdu_index_file(dl3Files,index_file_dir=dl3_dir)
    return

if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        config = yaml.safe_load(f)

    run_make_background(config)
