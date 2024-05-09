import os
import sys

import numpy as np
import yaml
from gammapy.data import DataStore


def prepare_dataset(config: dict, overwrite: bool = False) -> dict:
    """Prepare a dataset for analysis. Extract runs of interest from an existing datastore,
    reporting any missing runs.

    Parameters
    ----------
        config (dict)                              - dictionary with config information
        overwrite (bool)                           - bool to control whether or not to overwrite
                                                     the output directory
    Returns
    ----------
        config (dict)                              - dictionary with config information
                                                    (missing runs reported)

    """
    # Check if we have a list of runs of a file name
    runlist = config["run_selection"]["runlist"]
    if isinstance(runlist, str):
        with open(runlist, "r") as f:
            lines = f.readlines()
            runlist = np.array([int(line.strip()) for line in lines])

    # Open the data store
    db_dir = config["io"]["search_datastore"]
    data_store = DataStore.from_dir(db_dir)
    db_obs = data_store.obs_table["OBS_ID"]

    # Check which runs re in the data store
    obs_in_db = [ob for ob in runlist if ob in db_obs]
    missing_in_db = [ob for ob in runlist if ob not in db_obs]

    # Copy data of interest to a new dir
    in_dir = config["io"]["in_dir"]
    if not os.path.exists(in_dir):
        os.mkdir(in_dir)

    out_dir = config["io"]["out_dir"]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    results_dir = config["io"]["results_dir"]
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    # Copy
    try:
        data_store.copy_obs(obs_in_db, in_dir, overwrite=overwrite)
    except Exception as e:
        if len(obs_in_db) == 0:
            raise RuntimeError(
                f"Observations cannot be found in {db_dir}, do these files exist?"
            ) from e
        else:
            raise RuntimeError("Error copying files") from e

    # Update new config
    config["run_selection"]["runlist"] = obs_in_db
    config["run_selection"]["missing_runs"] = missing_in_db
    new_config = config["io"]["in_dir"] + "/config.yaml"
    with open(new_config, "w") as yamlfile:
        _ = yaml.dump(config, yamlfile)
        print(f"Written to {new_config}")

    return config


if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        config = yaml.safe_load(f)
