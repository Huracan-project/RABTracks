import huracanpy

# Calc
import numpy as np
import pandas as pd
import xarray as xr
from haversine import haversine
import geopandas as gpd

# Interface
from tqdm import tqdm
from glob import glob
import os
import pickle as pkl
import lzma

# Plots
import matplotlib.pyplot as plt
import seaborn as sns
import cartopy.crs as ccrs
from shapely.geometry import Polygon, MultiPolygon, Point
from matplotlib.lines import Line2D
from matplotlib.dates import HourLocator, DateFormatter, DayLocator
from matplotlib.patches import FancyBboxPatch
cm2inch = lambda L: [l*0.393701 for l in L]


# Analysis
from scipy.stats import linregress, pearsonr

# Parameters
BASINS = ["NA", "EP", "NI", "SA", "SI", "SP", "WP"]
SOURCES = ['IBTrACS', 
           'SyCLoPS-ERA5', 'SyCLoPS-JRA3Q', 'SyCLoPS-MERRA2',
       'TRACK-ECMWF-OP-AN', 'TRACK-ERA5', 'TRACK-JRA3Q', 'TRACK-MERRA2', 'TRACK-NCEP']

def clean(masked, dim = "record"):
    return masked.dropna(dim=dim, how="all")

# Plot parameters

cmap = sns.color_palette("tab20").as_hex()
PALETTE = {
    'IBTrACS':"black", 
    'TRACK-ERA5':cmap[0], # Dark blue
    'SyCLoPS-ERA5':cmap[1], # Light blue
    'TRACK-JRA3Q':cmap[2], # Dark orange
    'SyCLoPS-JRA3Q':cmap[3], # Light orange
    'TRACK-MERRA2':cmap[4], #Dark green
    'SyCLoPS-MERRA2':cmap[5], # Light green,
    'TRACK-ECMWF-OP-AN':cmap[6], # Violet
    'TRACK-NCEP':cmap[8], # Brown
}

ib_natures = ['DS', 'ET', 'MX', 'NR', 'SS', "TS",]
SyCLoPS_labels = ['DOTHL', 'DSD', 'DSE', 'DST', 'EX', 'HAL', 'HATHL', 'PL(PTLC)',
       'SC', 'SS(STLC)', 'TC', 'TD', 'TD(MD)', 'THL', 'TLO', 'TLO(ML)']

labels_colors = {
    #IBTrACS labels
    'DS':'lightyellow', 'ET':'dodgerblue', 'MX':'darkgrey',
    'NR':'lightgrey', 'SS':"green", "TS":"orange", 
    # SyCLoPS labels
    # Shallow LPSs or waves with weak surface circulations. DSD, DST and DSE are dry, tropical and extratropical DSs
    'DSD':"lightpink",
    'DSE':'lightblue',
    'DST':'khaki',
    'EX': 'dodgerblue', #ETC
    'HATHL': "darkviolet",#High altitude thermal low
    'SC': 'lightseagreen', # STC
    'SS(STLC)': "limegreen", # Subtropical Storm (Subtropical Tropical‐Like Cyclone)
    'TC':'red', # Tropical Cyclone
    'TD':'orange', # Tropical Depression
    'TLO':'gold' # Tropical Low
}
labels_hatches = {
    # IBTrACS labels
    'DS':'/', 'ET':'.', 'MX':'',
    'NR':'', 'SS':"/.", "TS":"//", 
    # SyCLoPS labels
    # Shallow LPSs or waves with weak surface circulations. DSD, DST and DSE are dry, tropical and extratropical DSs
    'DSD':"x",
    'DSE':'x',
    'DST':'x',
    'EX': '.', #ETC
    'HATHL': "+",#High altitude thermal low
    'SC': '/.', # STC
    'SS(STLC)': "/.", # Subtropical Storm (Subtropical Tropical‐Like Cyclone)
    'TC':'///', # Tropical Cyclone
    'TD':'//', # Tropical Depression
    'TLO':'/' # Tropical Low
}

ibtracs_units = {
    "sid": None,
    "season": "year",
    "number": None,
    "basin": None,
    "subbasin": None,
    "name": None,
    "iso_time": "UTC",
    "nature": None,
    "lat": "deg",
    "lon": "deg",
    "wmo_wind": "ms-1",
    "wmo_pres": "mb",
    "wmo_agency": None,
    "track_type": None,
    "dist2land": "km",
    "landfall": "km",
    "iflag": None,
    "usa_agency": None,
    "usa_atcf_id": None,
    "usa_lat": "deg",
    "usa_lon": "deg",
    "usa_record": None,
    "usa_status": None,
    "usa_wind": "ms-1",
    "usa_pres": "mb",
    "usa_sshs": None,
    "usa_r34_ne": "km",
    "usa_r34_se": "km",
    "usa_r34_sw": "km",
    "usa_r34_nw": "km",
    "usa_r50_ne": "km",
    "usa_r50_se": "km",
    "usa_r50_sw": "km",
    "usa_r50_nw": "km",
    "usa_r64_ne": "km",
    "usa_r64_se": "km",
    "usa_r64_sw": "km",
    "usa_r64_nw": "km",
    "usa_poci": "mb",
    "usa_roci": "km",
    "usa_rmw": "km",
    "usa_eye": "km",
    "tokyo_lat": "deg",
    "tokyo_lon": "deg",
    "tokyo_grade": None,
    "tokyo_wind": "ms-1",
    "tokyo_pres": "mb",
    "tokyo_r50_dir": None,
    "tokyo_r50_long": "km",
    "tokyo_r50_short": "km",
    "tokyo_r30_dir": None,
    "tokyo_r30_long": "km",
    "tokyo_r30_short": "km",
    "tokyo_land": None,
    "cma_lat": "deg",
    "cma_lon": "deg",
    "cma_cat": None,
    "cma_wind": "ms-1",
    "cma_pres": "mb",
    "hko_lat": "deg",
    "hko_lon": "deg",
    "hko_cat": None,
    "hko_wind": "ms-1",
    "hko_pres": "mb",
    "newdelhi_lat": "deg",
    "newdelhi_lon": "deg",
    "newdelhi_grade": None,
    "newdelhi_wind": "ms-1",
    "newdelhi_pres": "mb",
    "newdelhi_ci": None,
    "newdelhi_dp": "mb",
    "newdelhi_poci": "mb",
    "reunion_lat": "deg",
    "reunion_lon": "deg",
    "reunion_type": None,
    "reunion_wind": "ms-1",
    "reunion_pres": "mb",
    "reunion_tnum": None,
    "reuinion_ci": None,
    "reunion_rmw": "km",
    "reunion_r34_ne": "km",
    "reunion_r34_se": "km",
    "reunion_r34_sw": "km",
    "reunion_r34_nw": "km",
    "reunion_r50_ne": "km",
    "reunion_r50_se": "km",
    "reunion_r50_sw": "km",
    "reunion_r50_nw": "km",
    "reunion_r64_ne": "km",
    "reunion_r64_se": "km",
    "reunion_r64_sw": "km",
    "reunion_r64_nw": "km",
    "bom_lat": "deg",
    "bom_lon": "deg",
    "bom_type": None,
    "bom_wind": "ms-1",
    "bom_pres": "mb",
    "bom_tnum": None,
    "bom_ci": None,
    "bom_rmw": "km",
    "bom_r34_ne": "km",
    "bom_r34_se": "km",
    "bom_r34_sw": "km",
    "bom_r34_nw": "km",
    "bom_r50_ne": "km",
    "bom_r50_se": "km",
    "bom_r50_sw": "km",
    "bom_r50_nw": "km",
    "bom_r64_ne": "km",
    "bom_r64_se": "km",
    "bom_r64_sw": "km",
    "bom_r64_nw": "km",
    "bom_roci": "km",
    "bom_poci": "mb",
    "bom_eye": "km",
    "bom_pos_method": None,
    "bom_pres_method": None,
    "nadi_lat": "deg",
    "nadi_lon": "deg",
    "nadi_cat": None,
    "nadi_wind": "ms-1",
    "nadi_pres": "mb",
    "wellington_lat": "deg",
    "wellington_lon": "deg",
    "wellington_wind": "ms-1",
    "wellington_pres": "mb",
    "ds824_lat": "deg",
    "ds824_lon": "deg",
    "ds824_stage": None,
    "ds824_wind": "ms-1",
    "ds824_pres": "mb",
    "td9636_lat": "deg",
    "td9636_lon": "deg",
    "td9636_stage": None,
    "td9636_wind": "ms-1",
    "td9636_pres": "mb",
    "td9635_lat": "deg",
    "td9635_lon": "deg",
    "td9635_wind": "ms-1",
    "td9635_pres": "mb",
    "td9635_roci": "km",
    "neumann_lat": "deg",
    "neumann_lon": "deg",
    "neumann_class": None,
    "neumann_wind": "ms-1",
    "neumann_pres": "mb",
    "mlc_lat": "deg",
    "mlc_lon": "deg",
    "mlc_class": None,
    "mlc_wind": "ms-1",
    "mlc_pres": "mb",
    "usa_gust": "ms-1",
    "bom_gust": "ms-1",
    "bom_gust_per": "seconds",
    "reunion_gust": "ms-1",
    "reunion_gust_per": "seconds",
    "usa_seahgt": "m",
    "usa_searad_ne": "km",
    "usa_searad_se": "km",
    "usa_searad_sw": "km",
    "usa_searad_nw": "km",
    "storm_speed": "ms-1",
    "storm_dir": "degrees",
}