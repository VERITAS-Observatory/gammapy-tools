import os
from copy import deepcopy

import numpy as np

# Astropy stuff
from astropy import units as u
from astropy.table import Table
from gammapy.catalog import SourceCatalog3HWC, SourceCatalogGammaCat

# Gammapy stuff
from gammapy.irf import Background2D
from gammapy.maps import MapAxis
from gammapy.data import Observation
from scipy.ndimage import gaussian_filter


class BackgroundModelEstimator:
    """
    Largely based off https://docs.gammapy.org/0.18.2/tutorials/background_model.html


    """

    def __init__(
        self,
        energy: MapAxis,
        offset: MapAxis,
        smooth: bool = False,
        excluded_sources: list = [],
        smooth_sigma: float = 1,
    ) -> None:
        """BackgroundModelEstimator class for estimating the background from runs

        Parameters
        ----------
            energy (MapAxis)                    - Energy binning to be used
            offset (MapAxis)                    - Offset (to the camera centre) binning to be used
            smooth (bool)                       - Whether or not smooth the background
                                                 (default  False)
            excluded_sources (list)             - list of sources to be excluded
                                                 (astropy.coordinates.Skycoords)

        Returns
        ----------
            None

        """

        self.counts = self._make_bkg2d(energy, offset, unit="")
        self.exposure = self._make_bkg2d(energy, offset, unit="s TeV sr")
        self.excluded_sources = excluded_sources
        self.smooth = smooth
        self.smooth_sigma = smooth_sigma

        # Currelty hard coded...
        # ToDo wrap all this into package info
        this_dir, this_filename = os.path.split(__file__)
        try:
            star_path = os.path.join(this_dir, "Hipparcos_MAG8_1997.dat")
            self.star_data = np.loadtxt(star_path, usecols=(0, 1, 2, 3), skiprows=62)
        except Exception:
            star_path = os.path.join(
                os.environ.get("GAMMAPY_DATA"), "catalogs/", "Hipparcos_MAG8_1997.dat"
            )
            self.star_data = np.loadtxt(star_path, usecols=(0, 1, 2, 3), skiprows=62)
        self.star_cat = Table(
            {
                "ra": self.star_data[:, 0],
                "dec": self.star_data[:, 1],
                "id": self.star_data[:, 2],
                "mag": self.star_data[:, 3],
            }
        )

        # Assumes GAMMAPY_DATA is set
        self.cat = SourceCatalogGammaCat(
            "$GAMMAPY_DATA/catalogs/gammacat/gammacat.fits.gz"
        )
        self.hawc = SourceCatalog3HWC("$GAMMAPY_DATA/catalogs/3HWC.ecsv")

    @staticmethod
    def _make_bkg2d(energy: MapAxis, offset: MapAxis, unit: str) -> Background2D:
        """Get a 2D background with the desired axes

        Parameters
        ----------
            energy (MapAxis)                    - Energy binning
            offset (MapAxis)                    - Offset (to the camera centre) binning
            unit (str)                          - unit of the background

        Returns
        ----------
            bkg (Background2D)                  - Empty 2D background
        """
        return Background2D(axes=[energy, offset], unit=unit)

    def run(self, observations: list) -> None:
        """Generate background by stacking multiple backgrounds

        Parameters
        ----------
            observations (list)                 - Obserations for generating backgrounds


        Returns
        ----------
            None

        """
        for obs in observations:
            self.fill_counts(obs)

    #             self.fill_exposure(obs)

    def fill_counts(self, obs: Observation) -> None:
        """Fill the counts histograms for determining the background rate

        Parameters
        ----------
            obs (Observation)                   - Gammapy observation of events


        Returns
        ----------
            None

        """
        events = obs.events

        # Filter out regions of interest
        run_mask = self.exclude_known_sources(obs)
        run_mask = self.exclude_bright_stars(obs, run_mask=run_mask)

        energy_bins = self.counts.axes["energy"].edges
        offset_bins = self.counts.axes["offset"].edges

        # Bin the events
        counts = np.histogram2d(
            x=events.energy.to("TeV")[run_mask],
            y=events.offset.to("deg")[run_mask],
            bins=(energy_bins, offset_bins),
        )[0]

        self.counts.data += counts

        # keep this all one function
        # All events
        counts_all = np.histogram(events.offset.to("deg"), bins=offset_bins)[0] + 1e-9
        # Only kept events
        counts_exc = (
            np.histogram(events.offset[run_mask].to("deg"), bins=offset_bins)[0] + 1e-9
        )

        axes = self.exposure.axes
        offset = axes["offset"].center
        time = obs.observation_time_duration
        exposure = 2 * np.pi * offset * time * axes.bin_volume()

        # Scale exposure by fraction of events accepted
        self.exposure.quantity += exposure * (counts_exc / counts_all)

    # This could also be an exclusion file...
    def exclude_known_sources(
        self,
        obs: Observation,
        rad: float = 0.4,
        run_mask: np.array = None,
    ) -> np.array:
        """Exclude known sources from the background calculation

        Parameters
        ----------
            obs (Observation)                   - Gammapy observation of events
            rad (float)                         - radius of the exclusion region
            run_mask (np.array)                 - boolean array of if an event is to be included
                                                  (default None)


        Returns
        ----------
            run_mask (np.array)                 - boolean array of if an event is to be included

        """

        if run_mask is None:
            run_mask = np.ones(len(obs.events.radec.ra), dtype=bool)

        # Excluding Gammacat sources
        # Sources nearby
        gamma_cat_reduced_mask = (
            np.sqrt(
                (self.cat.table["ra"] - obs.fixed_pointing_info.radec.ra.deg) ** 2
                + (self.cat.table["dec"] - obs.fixed_pointing_info.radec.dec.deg) ** 2
            )
            < 2.5
        )

        # Apply Exclusion region to FoV
        for source in self.cat.table[gamma_cat_reduced_mask]:
            #             print (source["common_name"])
            # For extended source remove 3 times the extension?
            if source["morph_type"] in ["gauss", "shell"]:
                rad_ex = 3 * source["morph_sigma"] * u.deg
                # print (source["common_name"], rad_ex)
            else:
                rad_ex = 0.35 * u.deg

            run_mask *= (
                np.sqrt(
                    (obs.events.radec.ra - source["ra"] * u.deg) ** 2
                    + (obs.events.radec.dec - source["dec"] * u.deg) ** 2
                )
                > rad_ex
            )

        # Excluding HAWC sources
        # Sources nearby
        hawc_reduced_mask = (
            np.sqrt(
                (self.hawc.table["ra"] - obs.fixed_pointing_info.radec.ra.deg) ** 2
                + (self.hawc.table["dec"] - obs.fixed_pointing_info.radec.dec.deg) ** 2
            )
            < 2.5
        )

        # Apply Exclusion region to FoV
        # ToDo: Extended sources
        for source in self.hawc.table[hawc_reduced_mask]:
            run_mask *= (
                np.sqrt(
                    (obs.events.radec.ra - source["ra"] * u.deg) ** 2
                    + (obs.events.radec.dec - source["dec"] * u.deg) ** 2
                )
                > rad * u.deg
            )
        return run_mask

    # This could be sped up with a bright star file...
    def exclude_bright_stars(
        self,
        obs: Observation,
        rad: float = 0.35,
        mag: float = 8,
        run_mask: np.array = None,
    ) -> np.array:
        """Exclude bright stars from the background calculation

        Parameters
        ----------
            obs (Observation)                   - Gammapy observation of events
            rad (float)                         - radius of the exclusion region
                                                  default 0.35 deg
            mag (float)                         - magnitude below which stars are excluded
                                                  default 8.0
            run_mask (np.array)                 - boolean array of if an event is to be included
                                                  default None


        Returns
        ----------
            run_mask (np.array)                 - boolean array of if an event is to be included

        """
        if run_mask is None:
            run_mask = np.ones(len(obs.events.radec.ra), dtype=bool)

        # Look for stars above a mag cut and within the FoV
        srcs_mask = (
            np.sqrt(
                (self.star_cat["ra"] - obs.fixed_pointing_info.radec.ra.deg) ** 2
                + (self.star_cat["dec"] - obs.fixed_pointing_info.radec.dec.deg) ** 2
            )
            < 2.0
        )
        srcs_mask &= self.star_cat["mag"] < mag

        for src in self.star_cat[srcs_mask]:
            run_mask *= (
                np.sqrt(
                    (obs.events.radec.ra - src["ra"] * u.deg) ** 2
                    + (obs.events.radec.dec - src["dec"] * u.deg) ** 2
                )
                > rad * u.deg
            )
        return run_mask

    @property
    def background_rate(
        self,
    ) -> Background2D:
        """Return the background rate of the background

        Parameters
        ----------
            None


        Returns
        ----------
            rate (Background2D)                 - 2D background rate
        """
        rate = deepcopy(self.counts)

        rate.quantity /= self.exposure.quantity

        if self.smooth:
            rate = smooth(rate, sigma=self.smooth_sigma)

        return rate


def smooth(bkg: Background2D, sigma: float = 1.0) -> Background2D:
    """Smooths background rates from BackgroundModelEstimator.background_rate (bkg input)


    Parameters
    ----------
        bkg (Background2D)                      - 2D background to be smoothed
        sigma (float)                           - sigma of the gaussian for smoothing

    Returns
    ----------
        bkg (Background2D)                      - Smoothed 2D background

    """
    bkg_3d = bkg.to_3d()
    for i in range(len(bkg_3d.data)):
        smoothed = gaussian_filter(bkg_3d.data[i, :, :], sigma, 0)
        bkg_3d.data[i, :, :] = smoothed
    return bkg_3d.to_2d()
