# Background generation

Background generation is largely based off of [this outdated tutorial](https://docs.gammapy.org/0.18.2/tutorials/background_model.html).



## Generating background from multiple runs

When using config files, set `config["io"]["from_run"] = False` to generate a background by searching for runs within `config["io"]["search_datastore"]`. Runs matching the criteria set in `config["background_selection"]` will be included in the background calculation. Any known sources or bright stars will be excluded. Custom source exclusion isn't currently implemented (ToDo).


## Generating backgrounds from the source run

When using config files, set `config["io"]["from_run"] = True` to generate a background from the source run (the run we are interested). Using this option only the source run is used when generating backgrounds. If there are known sources and bright stars within the FoV they are excluded from the background calculation. Custom source exclusion isn't currently implemented (ToDo).



## Source Exclusion

By default, `BackgroundModelEstimator` will exclude events close to known VHE sources and bright stars. `BackgroundModelEstimator` queries GammaCat abnd 3HWC (as stored in Gammapy-data).

All sources with 2.5 degrees ofthe observation centre will be excluded. For point sources a circular exclusion region of radius 0.35 degrees is used. stars brighter than 8 mag have also have a 0.35 degree exclusion. For sources reported as extended in GammaCat or 3HWC, an exclusion region of 3 times the source extention is used (overkill)?

To correct the exposure for the excluded sources, the exposure is weighted by the ratio of all counts to the considered counts:
```python
        counts_all = np.histogram(
            events.offset.to("deg"),
            bins = offset_bins
        )[0] + 1e-9
        # Only kept events
        counts_exc = np.histogram(
            events.offset[run_mask].to("deg"),
            bins = offset_bins
        )[0] + 1e-9

        ...

        self.exposure.quantity += exposure * (counts_exc / counts_all)
```



## Smoothing

When using config files, set `config["background_selection"]["smooth"] = True` to "smooth" the background. This is done by converting the 2D background into a 3D background, smoothing it with a gaussian filter and then converting back to a 2D background.

```python
def smooth(bkg, sigma=1):
    """
    Smooths background rates from BackgroundModelEstimator.background_rate (bkg input)
    """
    bkg_3d = bkg.to_3d()
    for i in range(len(bkg_3d.data)):
        raw = bkg_3d.data[i,:,:]
        smoothed = gaussian_filter(bkg_3d.data[i,:,:],sigma,0)
        bkg_3d.data[i,:,:] = smoothed
    return bkg_3d.to_2d()
```