import sys
import yaml

# sys.path.insert(1,'../background_from_bkg/makeBackground')

# from BackgroundModelEstimator import *
# from makeBackground import *
# from makeRunwiseBackgrounds import *
from gammapy.data import DataStore
from ..make_background import get_background_for_run

config = sys.argv[1]
runwise = sys.argv[2]

with open(config, "r") as file:
    config = yaml.safe_load(file)

data_store = DataStore.from_dir(config["io"]["in_dir"])
obs_ids = data_store.obs_ids

for obs in obs_ids:
    print("Starting run # ", obs)
    get_background_for_run((obs, config))
