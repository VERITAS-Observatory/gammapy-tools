import logging
import numpy as np
from scipy.stats import norm
from os import environ
from IPython.display import display
import pandas as pd
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion
import matplotlib.pyplot as plt
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import MapDatasetOnOff, MapDataset,FluxPointsDataset
from gammapy.estimators import ExcessMapEstimator,FluxPointsEstimator, FluxPoints
from gammapy.makers import RingBackgroundMaker, SafeMaskMaker, MapDatasetMaker,SpectrumDatasetMaker
from gammapy.data import DataStore
from gammapy.maps import MapAxis, WcsGeom
from tqdm.auto import tqdm
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel
from astroquery.simbad import Simbad
from gammapy.maps import MapCoord
from gammapy.modeling import Fit
from astropy.table import Table
import sys
import seaborn as sns

sns.set_theme(font="Serif",font_scale=2,style='ticks',context='paper',palette='pastel')

datastore_dir = sys.argv[1]
source_pos = sys.argv[2] # [ra,dec] in deg 
plaw = sys.argv[3] # power law parameters as [index, amp] in TeV-1 s-1 cm-2
source = sys.argv[4]

logfile = f'{source}.log'

datastore = DataStore.from_dir(datastore_dir)
observations = datastore.get_observations()

source_pos = SkyCoord(source_pos[0],source_pos[1],unit='deg')

energy_axis = MapAxis.from_energy_edges([0.2,0.316,0.501,0.794,1.259,1.995,3.162,5.012,7.943,12.589,19.953,31.623]*u.TeV, unit="TeV")
energy_axis_true = MapAxis.from_energy_bounds(
    0.05, 100, 100, unit="TeV", name="energy_true"
)
geom = WcsGeom.create(
    skydir=(source_pos.ra.value, source_pos.dec.value),
    binsz=0.02,
    width=(4, 4),
    frame="icrs",
    proj="CAR",
    axes=[energy_axis],
)
simbad = Simbad()
simbad.reset_votable_fields()
simbad.add_votable_fields('ra', 'dec', "flux(B)", "flux(V)", "jp11")
simbad.remove_votable_fields('coordinates')

srcs_tab = simbad.query_region(source_pos, radius=1.5*u.deg)
srcs_tab = srcs_tab[srcs_tab["FLUX_B"]<6]
srcs_tab = srcs_tab[srcs_tab["FLUX_V"]!=np.ma.masked]

geom_image = geom.to_image().to_cube([energy_axis.squash()])

regions = CircleSkyRegion(center=source_pos, radius=0.5 * u.deg)
all_ex = [regions]
stars = []
for star in srcs_tab:
    pos = SkyCoord(star["RA"], star["DEC"], frame="fk5", unit=(u.hourangle, u.deg))
    star = CircleSkyRegion(center=pos, radius=0.3 * u.deg)
    stars.append(star)
    all_ex.append(star)

exclusion_mask = ~geom_image.region_mask(all_ex)
stacked = MapDataset.create(
    geom=geom, energy_axis_true=energy_axis_true, name="stacked"
)
offset_max = 1.75 * u.deg
maker = MapDatasetMaker()# add make exposure map?
maker.available_selection = []
maker_safe_mask = SafeMaskMaker(
    methods=["aeff-max","offset-max"], aeff_percent=10, offset_max = 1.75 * u.deg
)

for obs in observations:
    dataset = maker.run(stacked,obs)
    dataset = maker_safe_mask.run(dataset, obs)
    stacked.stack(dataset)

ring_maker = RingBackgroundMaker(
    r_in="0.63 deg", width="0.14 deg", exclusion_mask=exclusion_mask
)
dataset_on_off = ring_maker.run(stacked)

pl = PowerLawSpectralModel(index=plaw[0],amplitude=plaw[1]*u.Unit("TeV-1 s-1 cm-2"))
sky_model = SkyModel(spectral_model=pl,name='Crab')

on_region = CircleSkyRegion(center=source_pos, radius=np.sqrt(0.008) * u.deg)
spectral_dataset = dataset_on_off.to_spectrum_dataset(on_region,containment_correction=True)
spectral_dataset.models = [sky_model]
fit = Fit()
result = fit.run(datasets=spectral_dataset)

index = result.models['Crab'].index
index_err = result.models['Crab'].index.error
amplitude = result.models['Crab'].amplitude
amplitude_error = result.models['Crab'].amplitude

with open (logfile) as log:
    log.write(f'index: {index} +/- {index_err}')
    log.write(f'amplitude: {amplitude} +/- {amplitude_error}')

fpe = FluxPointsEstimator(
    energy_edges=[0.2,0.316,0.501,0.794,1.259,1.995,3.162,5.012,7.943,12.589]*u.TeV, source="Crab", selection_optional="all",reoptimize=False
)
flux_points = fpe.run(datasets=spectral_dataset)

flux_points.plot(sed_type='dnde',energy_power=2,color='plum')
result.models[0].spectral_model.plot([0.1,10]*u.TeV,sed_type='dnde',color='b',energy_power=2)
result.models[0].spectral_model.plot_error([0.1,10]*u.TeV,sed_type='dnde',alpha=0.1,energy_power=2)
plt.savefig(f'{source}_sed.pdf',format='pdf')

flux_points.write(f'{source}_fluxpoints.ecsv')

# make new dataset for significance calculations
stacked = MapDataset.create(
    geom=geom, energy_axis_true=energy_axis_true, name="stacked"
)

maker_safe_mask = SafeMaskMaker(
    methods=["offset-max"], offset_max = 1.75 * u.deg
)

for obs in observations:
    dataset = maker.run(stacked,obs)
    dataset = maker_safe_mask.run(dataset, obs)
    stacked.stack(dataset)

ring_maker = RingBackgroundMaker(
    r_in="0.63 deg", width="0.14 deg", exclusion_mask=exclusion_mask
)
dataset_on_off = ring_maker.run(stacked)

estimator = ExcessMapEstimator(np.sqrt(0.008) * u.deg, selection_optional=["alpha"],correlate_off=False,
                               #gamma_min_sensitivity=10,bkg_syst_fraction_sensitivity=0.05,
                               #apply_threshold_sensitivity=False,
                               #spectral_model = pl,
                               #sum_over_energy_groups=True,
                               )
lima_maps = estimator.run(dataset_on_off)

significance_map = lima_maps["sqrt_ts"]
excess_map = lima_maps["npred_excess"]

my_cmap = sns.color_palette("magma", as_cmap=True)

# We can plot the excess and significance maps
fig, (ax1, ax2) = plt.subplots(
    figsize=(13, 5), subplot_kw={"projection": lima_maps.geom.wcs}, ncols=2
)
ax1.set_title("Significance map")
significance_map.plot(ax=ax1, add_cbar=True,cmap=my_cmap,vmin=-5,vmax=5)
ax2.set_title("Excess map")
excess_map.plot(ax=ax2, add_cbar=True,cmap=my_cmap)
plt.savefig(f'{source}_sigmap.pdf', format='pdf')

significance_map_off = significance_map * exclusion_mask

# significance distribution
significance_all = significance_map.data[np.isfinite(significance_map.data)].flatten()
significance_off = significance_map_off.data[exclusion_mask & np.isfinite(significance_map_off.data)].flatten()

fig, ax = plt.subplots()
ax.hist(
    significance_all,
    density=True,
    alpha=0.5,
    color="red",
    label="all bins",
    bins=np.linspace(-5,10,50),
)

ax.hist(
    significance_off,
    density=True,
    alpha=0.5,
    color="blue",
    label="off bins",
    bins=np.linspace(-5,10,50),
)


mu, std = norm.fit(significance_off)
x = np.linspace(-5, 10, 100)
p = norm.pdf(x, mu, std)
ax.plot(x, p, lw=2, color="black")
ax.plot(x,norm.pdf(x,0,1),color='k',ls='--')

ax.legend()
ax.set_xlabel("Significance")
ax.set_yscale("log")
ax.set_ylim(1e-5, 1)
ax.set_xlim(-5, 10)

print(f"Fit results: mu = {mu:.2f}, std = {std:.2f}")
ax.text(-4.5, 0.5, f"Fit results: mu = {mu:.2f}, std = {std:.2f}")
plt.savefig(f'{source}_sigdist.pdf',format='pdf')

counts = lima_maps['counts'].get_by_coord((source_pos,1*u.TeV))[0]
sig = lima_maps['sqrt_ts'].get_by_coord((source_pos,1*u.TeV))[0]
bkg = lima_maps['npred_background'].get_by_coord((source_pos,1*u.TeV))[0]
alpha = lima_maps['alpha'].get_by_coord((source_pos,1*u.TeV))[0]

with open(logfile) as log:
    log.write(f'ON: {counts}')
    log.write(f'OFF: {bkg}')
    log.write(f'Significance: {sig}')
    log.write(f'Alpha: {alpha}')