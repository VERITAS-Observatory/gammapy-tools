import os
import numpy as np
import shutil
from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file
from astropy.coordinates import SkyCoord
from astropy import units as U
from gammapy.data import DataStore
from glob import glob

from .fake_location import LocationFaker
from ..make_background.background_tools import process_run


def find_runs(obs, config, n):

    # Get the KL divergence
    obs_table = process_run(
        obs,
        config,
        output_name=None,
        search_runs=None,
        bmimic=True,
        overwrite=False,
        njobs=config["config"]["njobs"],
    )

    # Sort by kl_div (reverse, kl should be small)
    obs_table.sort(["KL_DIV"], reverse=False)

    # Duration on the run of interest
    duration = obs_table["LIVETIME"][0]

    # The 0th item should be the same run
    obs_table = obs_table[1:]

    # Mask out runs within 10% livetime
    mask = (np.abs(obs_table["LIVETIME"] - duration) / duration) < 0.1

    # print (obs_table)
    obs_table = obs_table[mask]
    # print (obs_table)
    # Mask out listed bright sources:
    if "bright_sources" in config["background_selection"]:
        for source in config["background_selection"]["bright_sources"]:
            mask = obs_table["OBJECT"] == source
            obs_table = obs_table[~mask]

    # Check if we have enough observations remaining
    if len(obs_table) < n:
        n == len(obs_table)

    return obs_table[:n]


def mimic_data(config):

    runlist = config["run_selection"]["runlist"]
    # Default to 5 mimic datasets
    n_mimic = (
        config["background_selection"]["n_mimic"]
        if "n_mimic" in config["background_selection"]
        else 5
    )

    scrambe_theta = (
        config["background_selection"]["scramble_theta"]
        if "scramble_theta" in config["background_selection"]
        else 0.3
    )

    search_dir = config["io"]["search_datastore"]
    data_store = DataStore.from_dir(search_dir)

    # Run must have it's background so take data from out_dir
    input_dir = config["io"]["out_dir"]
    in_data = DataStore.from_dir(search_dir)

    # Make the output dirs
    for i in range(n_mimic):
        output_dir = input_dir + f"/mimic_{i + 1}"
        try:
            os.mkdir(output_dir)
        # if the directory already exists
        except OSError as error:
            print(error)

    faker = LocationFaker()
    for run in runlist:

        # Check if the run exists within the datastore
        of_interest = in_data.hdu_table["OBS_ID"] == run

        f_target = (
            search_dir
            + "/"
            + data_store.hdu_table[of_interest]["FILE_DIR"][0]
            + "/"
            + data_store.hdu_table[of_interest]["FILE_NAME"][0]
        )

        if not os.path.isfile(f_target):
            continue

        # print (n_mimic)
        mimic_runs = find_runs(run, config, n_mimic)
        # print (mimic_runs)
        # get random runs
        indx = np.random.randint(0, len(mimic_runs), n_mimic)
        for i in range(n_mimic):

            # Make sure file exists
            of_interest = (
                data_store.hdu_table["OBS_ID"] == mimic_runs["OBS_ID"][indx[i]]
            )

            f_mimic = (
                search_dir
                + "/"
                + data_store.hdu_table[of_interest]["FILE_DIR"][0]
                + "/"
                + data_store.hdu_table[of_interest]["FILE_NAME"][0]
            )
            if not os.path.isfile(f_mimic):
                continue

            # source_location
            target_location = SkyCoord(
                mimic_runs["RA_OBJ"][indx[i]] * U.deg,
                mimic_runs["DEC_OBJ"][indx[i]] * U.deg,
            )

            # Get the output name
            f_output = input_dir + f"/mimic_{i + 1}/" + os.path.basename(f_target)

            # todo add bright sources and stars
            known_sources = [target_location]

            faker.convert_fov(
                f_target,
                f_mimic,
                f_output,
                scramble_point=known_sources,
                overwrite=True,
                copy_background=True,
                scramble_theta=scrambe_theta,
            )

    # Make the datastores
    for i in range(n_mimic):

        out_dir = input_dir + f"/mimic_{i + 1}/"
        filelist = glob(out_dir + "/*.fits*")
        index_files = glob(out_dir + "/*index.fits*")
        filelist = list(set(filelist) - set(index_files))
        print(filelist)
        create_obs_hdu_index_file(filelist, out_dir)

        # Copy config
        shutil.copyfile(
            config["io"]["in_dir"] + "/config.yaml", out_dir + "/config.yaml"
        )
