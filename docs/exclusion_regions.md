# Excluding Regions

Regions can be automatically excluded using `gammapy_tools.utils.ExclusionFinder`. `ExclusionFinder` queries [gamma-cat](https://gamma-cat.readthedocs.io/), [3HWC](https://arxiv.org/abs/2007.08582) and the [Hipparcos star catalog](https://www.cosmos.esa.int/web/hipparcos/catalogues). `ExclusionFinder` queries local copies of the catalogs either installed via `gammapy-data` or included locally, therefore not requiring an external DB query.

## General Usage

The `ExclusionFinder` can find regions to exclude based on their distances to a given RA/Dec. For example:

```python
from gammapy_tools.utils import ExclusionFinder
from astropy.coordinates import SkyCoord

source_pos = SkyCoord.from_name("Crab")
...


exclusion = ExclusionFinder()
ex_regions, source_names = exclusion.find_exclusion(source_pos.ra.deg, source_pos.dec.deg, theta = 2., mag_cut=7)

```


In the above, we're excluding sources within 2 degrees of the Crab, with a cut on stars brighter than magnitude 7 (B band). This will find known VHE gamma-ray sources listed in either gamma-cat or 3HWC and bright stars.  This will return a list of (ra, dec, suggested exclusion radius), for example:

```
[(81.908688, 21.936965, 0.35),
 (84.411191, 21.142549, 0.35),
 (83.6279, 22.0243, 0.35),
 (85.166, 22.8719, 0.35),
 (83.63308, 22.0145, 0.35)]
```

If a source is listed as extended in gamma-cat, the suggested exclusion region with be 3 times the gaussian width of the source.

## Detailed Exclusion

`ExclusionFinder.find_exclusion` calls three functions under the hood which can also be accessed by the user.

`ExclusionFinder.find_gamma_sources` will query gamma-cat and return an `astropy.table.Table` with sources to be excluded and the following columns:
```
[
    "common_name", 
    "ra", 
    "dec", 
    "morph_type", 
    "morph_sigma"
]
```

Where `morph_type` is the morphological type and `morph_sigma` is the gaussian width if extended.


`ExclusionFinder.find_hawc_sources` will query 3HWC and return an `astropy.table.Table` with sources to be excluded and the following columns:
```
[
    "source_name",
    "ra", 
    "dec"
]
```

`ExclusionFinder.find_stars` will query Hipparcos star catalog and return an `astropy.table.Table` with sources to be excluded and the following columns:
```
[
    "id", 
    "mag", 
    "ra", 
    "dec"
]
```
Where `id` is the entry id in the table and `mag` is the B-band brightness of the star.

## API

### find_exclusion
::: gammapy_tools.utils.ExclusionFinder.find_exclusion
### find_gamma_sources
::: gammapy_tools.utils.ExclusionFinder.find_gamma_sources
### find_hawc_sources
::: gammapy_tools.utils.ExclusionFinder.find_hawc_sources
### find_stars
::: gammapy_tools.utils.ExclusionFinder.find_stars



