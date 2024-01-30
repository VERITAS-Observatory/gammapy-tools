# gammapy-tools
A repo for tools related to gammapy background generation.
This is the public version of [gammapy-background-analysis](https://github.com/VERITAS-Observatory/gammapy-background-analysis.git).


# Installation

This package requires [V2DL3](https://github.com/VERITAS-Observatory/V2DL3) and [Gammapy](https://gammapy.org/) to be installed.

Gammapy (v1.1) can be installed via:
```
pip install gammapy==1.1
```
Also install gammapy-data using
```
gammapy download datasets
```
And set the enviromental variable GAMMAPY_DATA, to that download location (for example):
```
export GAMMAPY_DATA=/path/to/gammapy/data
```

A modified version of `Hipparcos_MAG8_1997.dat` is used for finding stars. Copy this into the  `$GAMMAPY_DATA/catalogs/` directory:
```
cp background_from_bkg/Hipparcos_MAG8_1997.dat $GAMMAPY_DATA/catalogs/
```

As there maybe issues installing root-numpy, V2DL3 can be installed in one's enviroment as follows (replace `mamba` with `conda` for a conda install):
```
git clone https://github.com/VERITAS-Observatory/V2DL3.git
cd V2DL3

# Because its an env install of requirements
# Note the python version in the latest test throws issues with pytest
grep -A 100 "dependencies:" environment-eventdisplay.yml | grep "-" | grep -v "python" | awk '{print $2}' | xargs mamba install --yes

# root_numpy throws issues too. Only VEGAS uses it
mv setup.py _setup.py && grep -v "root_numpy" _setup.py > setup.py && pip install .
# clean up
cd ../ ; rm -r V2DL3
```

`background_from_bkg` can then be installed via:
```
pip install -r requirements.txt
pip install .
```

# Containerization
For details on working with containers see [containers.md](docs/containers.md)


# Usage

A number of worked examples can be found in [examples](examples).

To simplify analysis, a "driver" config file is used. An example of such a config file is shown in [config_crab.yaml](config_crab.yaml), which contains 4 crab runs.

* See [config.md](docs/config.md) for details on config files.
* See [mimicing_data.md](docs/mimicing_data.md) for details on mimicing and scrambling data.
