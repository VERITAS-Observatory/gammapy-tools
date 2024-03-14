import os
import numpy as np
from pyV2DL3.generateObsHduIndex import create_obs_hdu_index_file
from gammapy.data import DataStore
from glob import glob

from .fake_location import LocationFaker
from ..make_background.background_tools import process_run
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as U
import yaml


def find_runs(obs: int, config: dict, n: int) -> Table:
    """Find runs to use when mimicing the dataset

    Parameters
    ----------
        obs (int)                           - Observation ID for the observation of interest
        config (dict)                       - Configuration dictionary
        n (int)                             - Number of mimic runs to get


    Returns
    ----------
        obs_table (astropy.table.Table)     - Table of observations suitable for mimicking
    """

    # Get the KL divergence
    kl_file = config["io"]["out_dir"] + f"/{obs}_kl.fits"

    if os.path.isfile(kl_file):
        obs_table = Table.read(kl_file)
    else:
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
    # obs_table.sort(["KL_DIV"])

    # Duration on the run of interest
    duration = obs_table["LIVETIME"][0]

    # The 0th item should be the same run
    obs_table = obs_table[1:]

    # Mask out runs within 10% livetime
    mask = (np.abs(obs_table["LIVETIME"] - duration) / duration) < 0.1
    obs_table = obs_table[mask]

    # Mask out listed bright sources:
    if "bright_sources" in config["background_selection"]:
        for source in config["background_selection"]["bright_sources"]:
            mask = obs_table["OBJECT"] == source
            obs_table = obs_table[~mask]

    # Check if we have enough observations remaining
    if len(obs_table) < n:
        n == len(obs_table)

    return obs_table[:n]


def mimic_data(config: dict, randomise: bool = True) -> dict:
    """Creates mimic datasets for a runlist of interest

    Loops over each observation in the configuration file and finds the
    closest observations to be used when generating mimic datasets.
    Final datastore generation is randomized.

    Parameters
    ----------
        config (dict)                       - Configuration dictionary
        randomise (bool)                    - Whether to randomise the datasets
                                              if True (default) then a random suitable run is chosen
                                              if False then the most simiarly run
                                              is chosen for dataset 1, 2nd most for dataset 2,..,
                                              nth most for dataset n.

    Returns
    ----------
        config (dict)                       - Updated configuration dictionary
    """
    runlist = config["run_selection"]["runlist"]

    missing_runs = (
        config["run_selection"]["missing_runs"]
        if "missing_runs" in config["run_selection"]
        else []
    )

    # Remove missing runs from the runlist
    runlist = list(set(runlist) - set(missing_runs))

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

    questionable_runs = []
    for run in runlist:

        # Check if the run exists within the datastore
        of_interest = in_data.hdu_table["OBS_ID"] == run

        # f_target = (
        #     search_dir
        #     + "/"
        #     + data_store.hdu_table[of_interest]["FILE_DIR"][0]
        #     + "/"
        #     + data_store.hdu_table[of_interest]["FILE_NAME"][0]
        # )
        f_target = input_dir + os.path.basename(
            data_store.hdu_table[of_interest]["FILE_NAME"][0]
        )

        if not os.path.isfile(f_target):
            continue

        mimic_runs = find_runs(run, config, 10)

        if len(mimic_runs) < 1:
            questionable_runs.append((run, "no_mimic"))
            print(f"No mimic runs found for {run}")
            continue

        # Mask out large KL
        large_kl = mimic_runs["KL_DIV"] < 1.0
        mimic_runs = mimic_runs[large_kl]
        if len(mimic_runs) < 1:
            questionable_runs.append((run, "large_kl"))
            print(f"KL Value too large for {run}")
            continue

        indx = np.arange(len(mimic_runs))

        # Check if less than n_mimic runs are found
        if n_mimic > len(mimic_runs):
            n_repeat = n_mimic // len(mimic_runs)
            indx = np.repeat(indx, n_repeat + 1)

        # get random runs
        if randomise:
            # non repeating
            np.random.shuffle(indx)

        for i in range(n_mimic):

            # Make sure file exists
            of_interest = (
                data_store.hdu_table["OBS_ID"] == mimic_runs["OBS_ID"][indx[i]]
            )

            # print(
            #     f'Source Chosen: {mimic_runs["OBJECT"][indx[i]]}'
            #     + f' ({mimic_runs["OBS_ID"][indx[i]]}, kl = {mimic_runs["KL_DIV"][indx[i]]:0.3f})'
            # )

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
            # known_sources = []

            faker.convert_fov(
                f_target,
                f_mimic,
                f_output,
                scramble_point=known_sources,
                overwrite=True,
                copy_background=True,
                scramble_theta=scrambe_theta,
            )

    config["background_selection"]["questionable"] = questionable_runs

    # Make the datastores
    for i in range(n_mimic):

        out_dir = input_dir + f"/mimic_{i + 1}/"
        filelist = glob(out_dir + "/*.fits*")
        index_files = glob(out_dir + "/*index.fits*")
        filelist = list(set(filelist) - set(index_files))

        create_obs_hdu_index_file(filelist, out_dir)

        # Copy config
        # shutil.copyfile(
        #     config["io"]["in_dir"] + "/config.yaml", out_dir + "/config.yaml"
        # )
        with open(out_dir + "/config.yaml", "w") as outfile:
            yaml.dump(config, outfile, default_flow_style=False)

    if len(questionable_runs) > 0:

        print("Consider removing the following runs:")
        for i, (run, reason) in enumerate(questionable_runs):
            print(f"\t{run} ({reason})")

    return config
