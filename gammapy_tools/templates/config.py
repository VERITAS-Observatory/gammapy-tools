config = """io:
    # ED style runlist with one run number per line
    runlist: "/raid/romulus/swong/ED490/Eventdisplay_AnalysisScripts_VTS/scripts/crab_test.txt"
    # Where to search for background/mimic data
    search_datastore : "/raid/romulus/dl3/ed/dl3_fullenclosure_moderate2tel/"
    # Where to find your runlist's DL3 files (no bkg header)
    in_dir : "/raid/romulus/swong/mydl3/crabtest/"
    # Where to store your DL3 files with appended bkg - this is your DataStore for the analysis
    out_dir : "/raid/romulus/swong/mydl3/crabtestbkgs/"
    # Where to write results & plots
    results_dir: "/raid/romulus/swong/mydl3/crabtest/results/"

    # Whether or not to generate the background from the run
    from_run : false
    # Whether or not to analyze runs from a runlist or all runs in the results_dir
    from_runlist: true
    3d_bkg: false

binning:
    # Energy binning for bkg generation
    e_min: 0.1
    e_max: 100
    e_bins: 10
    # Offset/wobble binning for bkg generation
    off_min: 0
    off_max: 2.5
    off_bins: 5

config:
    # Number of cores to be used for parallel bkg generation
    njobs : 16

exclusion:
    # Star catalog for finding star exclusions
    star_cat : "./Hipparcos_MAG8_1997.dat"
    # Exculde region around a simbad resolveable source name
    # Expecting name, exclusion
    exclude_source_name :
    # Exculde region around a given ra, dec
    exclude_regions :
    # Exculde region around a given ra, dec
    # Exculde region around a given ra, dec
    # exclude_source_pos : [ ("287.10 deg", "6.39 deg", "0.5 deg"]
    # Exculde entire run within a distance to a given ra/dec
    # Expecting ra, dec, exclusion
    # exclude_source_run : [ ("287.10 deg", "6.39 deg", "3.5 deg"]

# Selection criteria for generating bkg headers
background_selection:
    # Maximum difference in NSB values
    nsb_diff : 1.5
    # Maximum 1/sin(el) difference
    el_diff : 1
    # Maximum azimuth difference
    az_diff : 45
    # Minimum number of telescopes
    n_tels : 4
    # Maximum time between runs in days - this doesn't need to be set if you're using same_szn
    time_max : 180
    # Whether bkgs should come from same season as data
    same_szn : true
    # Whether or not to smooth the background & by what width (in pixel space)
    smooth : true
    smooth_sigma : 1
    # How many hours of bkg headers to create
    time_request : 50
    # Whether to use KL Divergence criteria to pick nearest runs for bkgs
    # If time_request is set, this must be set to true
    KL_DIV: true
    # Whether or not to store the nearest bkgs as FITS files
    store_KL: false

run_selection:
    # Source of interest
    source_name : "Crab"
    # Whether to get source position from the DL3 header
    pos_from_DL3: true
    # If pos_from_DL3 is set to false, provide source coordinates in deg
    source_ra: 253.4675000
    source_dec: 39.7602778
    # Maximum offset in deg to select data from
    offset_max : 3.5
    # Minimum number of telescopes to select as valid runs
    n_tel : 3
    # Minimum elevation to select as valid runs
    el_min : 30

spectrum :
    # Energy threshold for spectral/integral flux calculations
    e_min : 0.2
    # Maximum energy for spectral/integral flux calculations
    e_max : 30
    # Number of energy bins for spectral calculations
    e_bins : 8
    # Spectral type
    type : "PL"
    # Spectral fit parameters (initial guess) - for power law: array of
    # [normalization, spectral index]
    # Assumes normalization is in units of TeV^-1 cm^-2 s^-1
    params : [3.45e-11,2.6]

lightcurve:
    #light curve bin size in minutes
    bin_size_min: 1000

sky_map:
    # Energy threshold for events to be included in the sky map (TeV)
    e_min: 0.01
    # Maximum energy for sky map (TeV)
    e_max: 100
    # Number of energy bins for the sky map
    n_bins: 30
    # Sky map size (deg)
    map_deg: 5
    # Sky map binning (deg) - ED default is 0.01, VEGAS is 0.02
    # Larger bins run faster
    bin_sz: 0.01
    offset_max: 1.75
    # Safe energy threshold (as a %) for effective area
    # - ED doesn't use this, so it's set fairly small
    aeff_max_percent: 0.1
    #exclusion regions are lists of lists of [[SkyCoord, radius], [SkyCoord, radius], ...]
    exclusion_regions: []
    #theta defines the ON region size (deg)
    theta: 0.0894427191
    #exclusion region around the source (deg)
    on_exclusion_region: 0.35
    #RBM parameters
    ring_rin: "0.6 deg"
    ring_width: "0.133333 deg"
    #whether or not to truncate the skymap colourbar at +/- 5 sigma
    truncate: true
    #dimmest visual magnitude for stars to be excluded
    min_star_brightness: 8

#where to save outputs
results_file: 'crab_test_results.yaml'
plot_names: 'crab_test'

"""
