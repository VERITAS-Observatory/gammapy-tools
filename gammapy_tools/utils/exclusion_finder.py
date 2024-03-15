from os import path, environ
from astropy.table import Table
import numpy as np
from gammapy.catalog import SourceCatalog3HWC, SourceCatalogGammaCat


class ExclusionFinder:

    def __init__(self):

        # ToDo wrap all this into package info
        this_dir, _ = path.split(__file__)
        try:
            star_path = path.join(this_dir, "../Hipparcos_MAG8_1997.dat")
            self.star_data = np.loadtxt(star_path, usecols=(0, 1, 2, 3), skiprows=62)

        except Exception:
            star_path = path.join(
                environ.get("GAMMAPY_DATA"), "catalogs/", "Hipparcos_MAG8_1997.dat"
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

    def find_gamma_sources(self, ra: float, dec: float, theta: float) -> Table:
        cat_mask = (
            (self.cat.table["ra"] - ra) ** 2 + (self.cat.table["dec"] - dec) ** 2
        ) < theta**2

        return self.cat.table[cat_mask][
            ["common_name", "ra", "dec", "morph_type", "morph_sigma"]
        ]

    def find_hawc_sources(self, ra: float, dec: float, theta: float) -> Table:
        hawc_mask = (
            (self.hawc.table["ra"] - ra) ** 2 + (self.hawc.table["dec"] - dec) ** 2
        ) < theta**2

        return self.hawc.table[hawc_mask][["source_name", "ra", "dec"]]

    def find_stars(
        self, ra: float, dec: float, theta: float, mag_cut: float = 7
    ) -> Table:
        star_mask = (
            (self.star_cat["ra"] - ra) ** 2 + (self.star_cat["dec"] - dec) ** 2
        ) < theta**2

        star_cat = self.star_cat[star_mask]

        return star_cat[star_cat["mag"] < mag_cut][["id", "mag", "ra", "dec"]]

    def find_exclusion(
        self,
        ra: float,
        dec: float,
        theta: float = 2.0,
        theta_cut: float = 0.35,
        mag_cut: float = 7,
    ):

        regions = []
        sources_excluded = []

        # Get the stars in the region
        stars = self.find_stars(ra, dec, theta, mag_cut)
        for star in stars:
            regions.append((star["ra"], star["dec"], theta_cut))
            sources_excluded.append(f"star_{star['id']:0.0f}")

        # Find the HAWC sources
        hawc = self.find_hawc_sources(ra, dec, theta)
        for source in hawc:
            regions.append((source["ra"], source["dec"], theta_cut))
            sources_excluded.append(source["source_name"])

        # Get the gammacat sources
        gamma = self.find_gamma_sources(ra, dec, theta)
        for source in gamma:
            theta_custom = theta_cut
            if source["morph_type"] in ["gauss", "shell"]:
                theta_custom = 3 * source["morph_sigma"]
            regions.append((source["ra"], source["dec"], theta_custom))
            sources_excluded.append(source["common_name"])

        return regions, sources_excluded

    def exclude_events(
        self,
        table: Table,
        ra: float,
        dec: float,
        theta: float = 2.0,
        theta_cut: float = 0.35,
        mag_cut: float = 7,
    ):
        exclude_regions, _ = self.find_exclusion(ra, dec, theta, theta_cut, mag_cut)
        mask = np.ones(len(table), dtype=bool)

        for reg in exclude_regions:
            mask &= ((table["RA"] - reg[0]) ** 2 + (table["DEC"] - reg[1]) ** 2) > reg[
                2
            ] ** 2

        return mask
