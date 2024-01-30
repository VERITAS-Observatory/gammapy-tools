import numpy as np
from astropy.time import Time


def time_cuts(observations, time_cut_file):
    time_intervals = []
    run, start, dur = np.genfromtxt(time_cut_file, unpack=True)
    for obs, i in enumerate(observations):
        if obs.id in run:
            start = Time(obs.gti.time_start[0], format="mjd")
            # unused?
            # SOB Added some place holder code
            stop = start + dur
            time_intervals.append([start, stop])
