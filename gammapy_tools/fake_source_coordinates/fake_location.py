import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d

from ..utils import get_cdf

# from copy import deepcopy


class LocationFaker:
    """Class to change the coordiates within a FoV
    Change camara X/Y to a different RA/Dec for a matched set of observations

    """

    def __init__(self) -> None:
        pass

    def load_target(self, filename: str) -> None:
        """Load the target FoV from a DL3 file

        Parameters
        ----------
            filename (str)                      -  dl3 file name

        Returns
        ----------
            None
        """

        self.target_file = filename
        with fits.open(filename) as target_hdul:
            self.target_ra = target_hdul["EVENTS"].header["RA_PNT"]
            self.target_dec = target_hdul["EVENTS"].header["DEC_PNT"]
            self.target_obj_ra = target_hdul["EVENTS"].header["RA_OBJ"]
            self.target_obj_dec = target_hdul["EVENTS"].header["DEC_OBJ"]
            self.obs_id = target_hdul["EVENTS"].header["OBS_ID"]

    def load_faked(self, filename: str) -> None:
        """Load the faked FoV from a DL3 file

        Parameters
        ----------
            filename (str)                      - dl3 file name

        Returns
        ----------
            None
        """

        self.faked_file = filename
        with fits.open(filename) as faked_hdul:
            self.faked_ra = faked_hdul["EVENTS"].header["RA_PNT"]
            self.faked_dec = faked_hdul["EVENTS"].header["DEC_PNT"]
            self.faked_obj_ra = faked_hdul["EVENTS"].header["RA_OBJ"]
            self.faked_obj_dec = faked_hdul["EVENTS"].header["DEC_OBJ"]
            self.faked_coords = SkyCoord(
                self.faked_ra, self.faked_dec, unit="deg", frame="icrs"
            )

    def convert_fov(
        self,
        source_file: str,
        faked_file: str,
        output_name: str,
        overwrite: bool = False,
        copy_background: bool = False,
        scramble_point: list = None,
        scramble_theta: float = 0.3,
    ) -> None:
        """Convert the coordinates of the FoV

        Parameters
        ----------
            source_file (str)                   - dl3 file name for the source fov
            faked_file (str)                    - dl3 file name for the faked fov
            output_name (str)                   - dl3 file name for the output file
            overwrite (bool)                    - whether or not to overwrite output_name
                                                  (default False, do not overwrite)
            copy_background (bool)              - whether or not to copy the background
                                                  from source_file
                                                  (default False, output_name contains
                                                   faked_file's background)
            scramble_point (list)               - Points to be scrambled list of astropy
                                                  sky coordinates
                                                  (default None)
            scramble_theta (float or list)      - Radius of the region to be scrambled and
                                                  half width of anulus region of scrambling
                                                  (r_outter - r_inner = 2 * scramble_theta)



        Returns
        ----------
            None
        """

        # Load in the files
        self.load_target(source_file)
        self.load_faked(faked_file)

        # Determine RA/Dec offset
        ra_offset = self.target_ra - self.faked_ra
        dec_offset = self.target_dec - self.faked_dec

        faked_hdul = fits.open(faked_file)
        faked_table = Table.read(faked_hdul["EVENTS"])

        if scramble_point is not None:
            x_to_ra = self.faked_ra
            y_to_dec = self.faked_dec
            faked_table = self.scramble_data(
                faked_table, scramble_point, scramble_theta, x_to_ra, y_to_dec
            )

        faked_table["RA"] += ra_offset
        faked_table["DEC"] += dec_offset

        # Reassign the modifed data
        faked_hdul["EVENTS"].data = fits.BinTableHDU(faked_table).data

        # Modify the header
        faked_hdul["EVENTS"].header["RA_PNT"] = self.target_ra
        faked_hdul["EVENTS"].header["DEC_PNT"] = self.target_dec

        # Add in the original ra/dec
        faked_hdul["EVENTS"].header["RA_ORG"] = (
            self.faked_ra,
            "RA pointing of the original run",
        )
        faked_hdul["EVENTS"].header["DEC_ORG"] = (
            self.faked_dec,
            "DEC pointing of the original run",
        )
        # Add in the original ra/dec of the source
        faked_hdul["EVENTS"].header["RA_SRC"] = (
            self.faked_obj_ra,
            "RA of the original object",
        )
        faked_hdul["EVENTS"].header["DEC_SRC"] = (
            self.faked_obj_dec,
            "DEC of the original object",
        )
        # Add in the ra/dec offset
        faked_hdul["EVENTS"].header["RA_OFF"] = (ra_offset, "RA offset from original")
        faked_hdul["EVENTS"].header["DEC_OFF"] = (
            dec_offset,
            "DEC offset from original",
        )

        # Change the obs_id
        faked_hdul["EVENTS"].header["OBS_ID"] = self.obs_id
        faked_hdul["EFFECTIVE AREA"].header["OBS_ID"] = self.obs_id

        # Copy background
        target_hdul = fits.open(source_file)

        if copy_background:
            try:

                if "BACKGROUND" in faked_hdul:
                    bkg_indx = faked_hdul.index_of("BACKGROUND")
                    faked_hdul.pop(bkg_indx)
                # if "BACKGROUND" in faked_hdul:
                #     faked_hdul["BACKGROUND"] = target_hdul["BACKGROUND"]
                # else:
                faked_hdul.append(target_hdul["BACKGROUND"])
                faked_hdul[-1].name = "BACKGROUND"

            except Exception as e:
                print("Issue writing background")
                print(e)

        faked_hdul.writeto(output_name, overwrite=overwrite)
        target_hdul.close()

    def scramble_data(
        self,
        tab: Table,
        scramble_source: list,
        scramble_theta: float,
        ra_off: float,
        dec_off: float,
    ) -> Table:
        """Scramble events in a ring around the centre based on the distance to the source

        Parameters
        ----------
            tab (astropy.table.Table)           - Table of DL3 EVENTS
            scramble_source (list)              - list of sources to scramble about
            scramble_theta (float or list)      - radius of region around source to be scrambled
            ra_off (float)                      - Right Ascension offset (deg)
            dec_off (float)                     - Declination offset (deg)

        Returns
        ----------
            tab (astropy.table.Table)           - Table of scrambled DL3 EVENTS
        """

        # if only a single theta is given, cast it to a list
        if isinstance(scramble_theta, float):
            scramble_theta = [scramble_theta] * len(scramble_source)
        else:
            if len(scramble_theta) != len(scramble_source):

                err_str = (
                    "Expecting scramble_source and scramble theta to be of equal length"
                )
                err_str += f"\n({len(scramble_source)}, {len(scramble_theta)})"
                raise TypeError(err_str)

        dists = []
        rnd_angs = []
        rng_dists = []
        # Loop over sources
        for source, theta in zip(scramble_source, scramble_theta):
            # Calculating the distance to the source
            dist = self.faked_coords.separation(source).to("deg").value
            dists.append(dist)

            # Mask of events within the radial distance to the source of interest
            event_mask = (
                np.abs(tab["Xoff"] ** 2 + tab["Yoff"] ** 2 - dist**2) < theta**2
            )

            # Masking events around the source of interest
            src_mask = (
                np.abs(
                    (tab["RA"] - self.faked_obj_ra) ** 2
                    + (tab["DEC"] - self.faked_obj_dec) ** 2
                )
                < theta**2
            )

            # Get a random angle [radians]
            rnd_ang = np.random.uniform(
                low=0, high=2 * np.pi, size=len(tab[event_mask])
            )
            rnd_angs.append(rnd_ang)

            # Random distance between +/- theta/2
            # This shouldn't be uniform, weight by the radial acceptances across +/- theta/2
            # Sample the radial distribtuion out of the source region
            bw = 0.05
            rad_bins = np.arange(
                dist - theta / 2.0 - bw, dist + theta / 2.0 + 1.1 * bw, bw
            )
            rad_binsc = rad_bins[:-1] + 0.5 * (rad_bins[1:] - rad_bins[:-1])

            # Get radial distribution from outside of the source region
            radial, _ = np.histogram(
                np.sqrt(tab["Xoff"][~src_mask] ** 2 + tab["Yoff"][~src_mask] ** 2),
                bins=rad_bins,
            )
            radial_mask = (rad_binsc > (dist - theta / 2.0)) & (
                rad_binsc < (dist + theta / 2.0)
            )

            rad_cdf = get_cdf(radial[radial_mask])
            # Use an interpolation to get any probability
            inter_rad = interp1d(
                rad_cdf,
                rad_binsc[radial_mask],
                fill_value="extrapolate",
                bounds_error=False,
            )
            rnd_num = np.random.uniform(0, 1, size=len(tab[event_mask]))

            rnd_dist = inter_rad(rnd_num)
            rng_dists.append(rnd_dist)

            # Getting random x/y
            rnd_x = (rnd_dist) * np.cos(rnd_ang)
            rnd_y = (rnd_dist) * np.sin(rnd_ang)

            # Reassign the camera coordinates
            tab["Xoff"][event_mask] = rnd_x
            tab["Yoff"][event_mask] = rnd_y

            # Reassign the sky coordinates
            tab["RA"][event_mask] = tab["Xoff"][event_mask] + ra_off
            tab["DEC"][event_mask] = tab["Yoff"][event_mask] + dec_off

        return tab
