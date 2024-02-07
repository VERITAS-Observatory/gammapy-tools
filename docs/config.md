# Config Files

## io Section

The `io` section defines where files are stored.
```
io:
    search_datastore : "/gamma_data/MegaStore/moderate2tel/"
    in_dir : "/gamma_data/Crab_data/"
    out_dir : "/gamma_data/Crab_data_background/"
    results_dir : "/gamma_data/Crab_data_background/results/"
    # Whether or not to generate the background from the run
    from_run : False  
```

* `search_datastore` is the path to a datastore that is queried when trying to find suitable background runs.

* `in_dir` is the input directory/datastore, where the unmodified fits files will be located. Backgrounds will be obtained for run in this directory.

* `out_dir` is where runs with background files will be saved to. `out_dir` != `in_dir` to avoid overwriting raw data. However since backgrounds can overwritten without effecting raw event data, these could be the same directory.

* `results_dir` is where high level results (light curves, spectra, skymaps, etc.) will be written to.

* `from_run` is a boolean (True/False) to control whether or not backgrounds are generated from the run. 

    * If `from_run == True` then backgrounds are generated from the FoV of each run, with known sources and bright stars excluded. Backgrounds will only be generated from 1 run.
    
    * If  `from_run == False` then `search_datastore` will be queried for suitable background runs. Known sources and bright stars will be excluded. Background will be generated from `n` runs matching the selection criteria.

## binning Section

The `binning` section defines the bins used when generating backgrounds. 

```
binning:
    # Energy binning
    e_min: 0.1
    e_max: 100
    e_bins: 20
    # Offset/wobble binning
    off_min: 0
    off_max: 2.5
    off_bins: 5
```

* `e_min` the minimum energy (TeV) 
* `e_max` the maximum energy (TeV) 
* `e_bins` the number of bins between `e_min` and `e_max`
* `off_min` the minimum offset to the camera centre (degrees)
* `off_max` the maximum offset to the camera centre (degrees)
* `off_bins` the number of bins between `off_min` and `off_max`



## config Section

The `config` section defines some run parameters relevent to the execution of scripts.

```
config:
    njobs : 10
```

* `njobs` the number of parallel jobs allowed. 
    * `njobs == -1` will use all available cores
    * `njobs == 1` will use a single core



## exclusion Section
The `exclusion` section defines regions to be excluded when generating backgrounds and running automatic analysis.

```
exclusion:
    star_cat : "./Hipparcos_MAG8_1997.dat"
    # Exculde region around a given ra, dec
    exclude_regions :
        - ["287.10", "6.39", "0.5"]
    # Exculde region around a simbad resolveable source name
    # Expecting name, exclusion
    exclude_source_name : [ "MGRO J1908+06", "0.5 deg"]
```


* `star_cat` path to the star catalog used (redundent?)
* `exclude_regions` list of regions to exclude. Each entry is assumed to be `[ra (deg), dec (deg), exclusion_radius (deg)]`. Where `ra` and `dec` are the Right Ascension and Declination in degrees (J2000). `exclusion_radius` is the radius of a circular region to exclude. 



## background_selection Section

The `background_selection` defines the search parameter to be used when searching for background runs.

```
background_selection:
    # el_min : 40
    # el_max : 30
    nsb_diff : 0.7
    # 1/sin(el) difference
    el_diff : 0.15
    # az difference
    az_diff : 45
    # Minimum number of telescopes
    n_tels : 4
    # Maximum time between runs in days
    time_max : 180
    # Whether or not to smooth the background
    smooth : False
    # Runlist for creating backgroud
    # bkg_runlist : 
    #     64080 : [65895, 63978]
    #     64081 : [65895, 63978]
    #     64082 : [65895, 63978]
    #     64083 : [65895, 63978]
```

* `el_min` minimum elevation for background runs.
* `el_max` maximum elevation for background runs.
* `el_dif` maximum different in 1/sin(elevation).
* `nsb_dif` maximum difference between nsb levels.
* `az_dif` maximum difference in azimutal angle.
* `n_tels` minimum number of telescopes
* `time_max` maximum time difference between runs (days)
* `smooth` whether or not to smooth the background:
    * If `smooth == True` 2D background is converted to 3D, smoothed with a gaussian filter and converted back to 2D.
    * If `smooth == False` the 2D background isn't smoothed.
* `bkg_runlist` (optional) if `bkg_runlist` is provided, then the runs explicitly listed will be used when generating backgrounds. (Not currently implemented... ToDo: check why not)


# analysis_selection Section

The `analysis_selection` section defines analysis run parameters for running automatic analysis. (ToDo: Add energy, offset cuts and any other analysis parameters).

```
analysis_selection:
    theta2 : 0.008
```

* `theta2` is the $\theta^2$-cut used when running analysis.


# run_selection Section

The `run_selection` select defines the search parameters when searching for data we're interested in.

```
run_selection:
    # Runs of interests
    runlist : [ 7, 64080, 64081, 64082, 64083]    
    # Source of interest
    source_name : "Crab"
    source_ra : 83.6333
    source_dec : 22.0133
```

* `runlist` array of runs or a text file runlist
* `source_name` name of the source
* `source_ra` Right Ascension in degrees (J2000)
* `source_dec`  Declination in degrees (J2000)


`runlist` can contain runs that are not included in `io:search_datastore` (for example in the above `7`). When running `prepare_dataset`, missing runs will be reported as `run_selection:missing_runs`.

Automatic searching based on source name or coordinates isn't yet...