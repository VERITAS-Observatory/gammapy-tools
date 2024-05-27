## 1.0.1 (2024-05-09)

### Fix

- **pyproject.toml**: fixing versioning
- **gammapy_tools/make_background/prepare_data.py**: adding overwrite option

### Refactor

- **gammapy_tools/__version__.py**: post test update
- **__version__**: switching to __version__ for the version number

## 1.0.0 (2024-05-07)

### Feat

- **utils/exclusion_finder.py**: adding exclusion finder
- **gammapy_tools/fake_source_coordinates/process.py**: adding mimic_data
- **make_background**: parallel reduction of background and mimic search
- **gammapy_tools/make_background**: implementing a parallel + reduction method to spead up background generation
- **make_background**: implementation of closed N background runs
- **background_models**: adding a user defined smoothing sigma

### Fix

- **pyproject.toml**: changing gammapy version
- **gammapy_tools/make_background/background_models.py**: updating to gammapy1.2
- **gammapy_tools/utils/exclusion_finder.py**: adding a check for gammacat and hawc
- **gammapy_tools/utils/exclusion_finder.py**: correcting path
- **gammapy_tools/utils/exclusion_finder.py**: correcting imports and file paths
- **fake_source_coordinates/process.py**: adding safety check
- **gammapy_tools/fake_source_coordinates/process.py**: fixing print and removing target
- **process.py**: correcting background
- **make_background**: storing kl_div table
- **background_tools.py**: exposure searching
- **prepare_data.py**: adding more info
- **prepare_data.py**: raise error when no runs are found
- **gammapy_tools/analysis/rbm.py**: removing hard coded map size values
- **Hipparcos_MAG8_1997.dat**: updated hipparcos file to have correct colour column
- **analysis**: update the analysis scripts to work with config files

### Refactor

- **gammapy_tools/fake_source_coordinates/process.py**: adding more robust file finder
- **fake_source_coordinates/process.py**: reincluding source when scrambling
- **gammapy_tools/fake_source_coordinates/fake_location.py**: popping exisiting background
- **make_background**: removing ignore warnings
- **make_background.py**: removing additional background and allow previous method
- **templates/config.py**: adding KL_DIV
- **utils**: moving functions to utils
- **analysis_notebook**: updating example
- **data_products**: removing plt.show
- **_version.**: bumping version
- **analysis_notebook**: migrating analysis notebook
- **fake_source_coordinates**: adding __all__
- **pyproject.toml**: changing required to minimum python version

### Perf

- **gammapy_tools/utils/run_details.py**: Changed the mimic selection to perform a nested search

## 0.1.2 (2024-01-30)

### Fix

- **gammapy_tools/analysis**: Adding initial analysis scripts
- **templates**: get_config
- **templates**: Adding an analysis config template

### Refactor

- **repository**: running pre-commit linting and formating on previous commit
- **background_from_bkg**: adding initial code
