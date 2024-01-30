config = """io:
    search_datastore : "/path/to/datastore/"
    in_dir : "/path/to/input_dir/"
    out_dir : "/path/to/output_dir/"
    results_dir : "/path/to/results_dir/"
    # Whether or not to generate the background from the run
    from_run : False

binning:
    # Energy binning
    e_min: 0.1
    e_max: 100
    e_bins: 20
    # Offset/wobble binning
    off_min: 0
    off_max: 2.5
    off_bins: 5

config:
    njobs : 1

exclusion:
    star_cat : "./Hipparcos_MAG8_1997.dat"
    # Exculde region around a simbad resolveable source name
    # Expecting name, exclusion
    exclude_source_name :
        - [ "MGRO J1908+06", "0.5 deg"]
    # Exculde region around a given ra, dec
    exclude_regions :
        - ["287.10", "6.39", "0.5"]
    # Exculde region around a given ra, dec
    # Exculde region around a given ra, dec
    # exclude_source_pos : [ ("287.10 deg", "6.39 deg", "0.5 deg"]
    # Exculde entire run within a distance to a given ra/dec
    # Expecting ra, dec, exclusion
    # exclude_source_run : [ ("287.10 deg", "6.39 deg", "3.5 deg"]


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

analysis_selection:
    theta2 : 0.008

run_selection:
    # Runs of interests
    runlist : [ 64080, 64081, 64082, 64083]
    # Source of interest
    source_name : "Crab"
    source_ra : 83.6333
    source_dec : 22.0133
    # Coords to search around
    # coords : ["287.10 deg", "6.39 deg"]
    # coords_galactic : ["40.39 deg", "-0.79 deg"]
    # Maximum offset to select data from
    offset_max : 1.5
    # Number of telescopes to select
    n_tel : 4
    # elevation min/maz
    el_min : 40
    # el_max : 90
    # Time range (d/m/y)
    # time_min : "01/09/2007"
    # time_max : "01/09/2025"
"""