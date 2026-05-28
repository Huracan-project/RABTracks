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

def save_xarray_to_pickle(obj, path, compress = True):
    """
    Save an xarray Dataset/DataArray to a pickle file in a way that's portable:
    - Materialize all data into memory (no lazy netCDF/dask-backed arrays)
    - Detach from any file handles
    - Remove filesystem path hints like encoding["source"]
    - Verify materialization by forcing numpy conversion of all variables

    Parameters
    ----------
    obj : xarray.Dataset or xarray.DataArray
        The xarray object to save.
    path : str or pathlib.Path
        Output pickle path.
    """
    if not isinstance(obj, (xr.Dataset, xr.DataArray)):
        raise TypeError(f"obj must be an xarray Dataset or DataArray, got {type(obj)!r}")

    # Work on a detached copy (avoid mutating caller's object)
    # deep=True ensures we own the arrays once loaded.
    obj2 = obj.copy(deep=True)

    # Materialize data (handles both netCDF-backed and (most) dask-backed cases)
    # - .load() is xarray-native and will compute dask if present.
    obj2 = obj2.load()

    # Close any open file handles (safe no-op if nothing open)
    try:
        obj2.close()
    except Exception:
        pass

    # Strip path-ish encodings that can leak local filesystem info
    def _strip_encodings(ds):
        if hasattr(ds, "encoding") and isinstance(ds.encoding, dict):
            ds.encoding.pop("source", None)
            ds.encoding.pop("original_shape", None)

        # Dataset has .variables; DataArray does not in the same way
        if isinstance(ds, xr.Dataset):
            for v in ds.variables:
                enc = getattr(ds[v], "encoding", None)
                if isinstance(enc, dict):
                    enc.pop("source", None)
                    enc.pop("original_shape", None)

    if isinstance(obj2, xr.Dataset):
        _strip_encodings(obj2)
    else:  # DataArray
        _strip_encodings(obj2)
        # Also strip the underlying Variable encoding if present
        try:
            if isinstance(obj2.variable.encoding, dict):
                obj2.variable.encoding.pop("source", None)
                obj2.variable.encoding.pop("original_shape", None)
        except Exception:
            pass

    # Verify: force actual data access; if anything is still lazy, this is where it breaks.
    if isinstance(obj2, xr.Dataset):
        for name, var in obj2.variables.items():
            _ = np.asarray(var.data)
        # Stronger forcing option (can be expensive); uncomment if desired:
        # _ = obj2.to_dataframe()
    else:
        _ = np.asarray(obj2.data)

    # Write pickle
    with open(path, "wb") as f:
        pkl.dump(obj2, f,)

    if compress:
        with lzma.open(path+".xz", "wb") as f:
            pkl.dump(obj2, f)
    

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
    'HAL': "violet", # High‐altitude Low
    'HATHL': "darkviolet",#High altitude thermal low
    'SC': 'lightseagreen', # STC
    'SS(STLC)': "limegreen", # Subtropical Storm (Subtropical Tropical‐Like Cyclone)
    'PL(PTLC)': "turquoise", # Polar Low (Polar Tropical‐Like Cyclone)
    'HTC':'yellowgreen', # Hybrid TC
    'TC':'orange', # Tropical Cyclone
    'TD':'gold', # Tropical Depression
    'TLO':'lemonchiffon' # Tropical Low
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
    'HAL': "+", # High‐altitude Low
    'HATHL': "++",#High altitude thermal low
    'SC': '/..', # STC
    'SS(STLC)': "/.", # Subtropical Storm (Subtropical Tropical‐Like Cyclone)
    'PL(PTLC)': "//..", # Polar Low (Polar Tropical‐Like Cyclone)
    'HTC':"//.", # Hybrid TCs
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

from numpy.typing import ArrayLike
from itertools import groupby
# Function from Leo Saffin https://github.com/Huracan-project/wcsi/blob/main/wcsi/nature.py#L252
def nature(
    b: ArrayLike,
    vtl: ArrayLike,
    vtu: ArrayLike,
    vort: ArrayLike,
    is_tc: ArrayLike,
    *,
    b_threshold: float = 15,
    vtl_threshold: float = 0,
    vtu_threshold: float = 0,
    vort_threshold: float = 6,
    min_count: int = 4,
    et: bool = False,
    smooth: bool = False,
) -> np.ndarray:
    """Derive a nature tag from the cyclone structure

    TC - Tropical Cyclone
    Vo - Weak Vortex
    BC - Baroclinic
    Tr - Trough
    MV - Mid-level vortex
    Ot - Other

    If et=True
    ET - Extratropical transition
    WS - Warm seclusion

    Parameters
    ----------
    b
        Cyclone phase space asymmetry
    vtl
        Cyclone phase space low-level warm core
    vtu
        Cyclone phase space upper-level warm core
    vort
        850hPa vorticity
    is_tc
        Points previously used to identify the cyclone as tropical cyclone (e.g. WCSI)
    b_threshold
        The threshold of the asymmetry parameter, below which is considered to be a
        tropical cyclone
    vtl_threshold
        The threshold of the low-level warm-core parameter, above which is considered
        to be a tropical cyclone
    vtu_threshold
        The threshold of the upper-level warm-core parameter, above which is considered
        to be a tropical cyclone
    vort_threshold
        The minimum threshold for 850-hPa vorticity, below which is considered as a
        weak vortex
    min_count
        Number of
    et
        Add labels for extratropical transition. Each stage must last for min_count
        points
    smooth
        Add a smoothing to the nature tags. Any excursions less than min_count are
        removed

    Returns
    -------


    """
    nat = np.zeros(len(vort), dtype="U2")

    # Too weak = vortex
    weak = vort < vort_threshold
    nat[weak] = "Vo"

    # WCSI label as tropical cyclone
    nat[is_tc] = "TC"

    # Other CPS categories
    symmetric = b <= b_threshold
    warm_core = vtl > vtl_threshold
    trough = vtu <= vtu_threshold

    # Any Warm core/symmetric periods adjacent to TC are also TC
    # Label as tropical storm for now
    nat[(nat == "") & (b <= b_threshold) & (vtl > vtl_threshold)] = "TS"
    nat_consecutive = [(k, sum(1 for _ in g)) for k, g in groupby(nat)]
    idx = 0
    for m, (nat_, count) in enumerate(nat_consecutive):
        if nat_ == "TS":
            # Allow for <1 day excursions between TC-TS
            idx_start = max(0, idx - min_count)
            if (nat[idx_start : idx + count] == "TC").any():
                nat[idx_start : idx + count] = "TC"

            idx_end = min(len(nat), idx + count + min_count)
            if (nat[idx + count : idx_end] == "TC").any():
                nat[idx:idx_end] = "TC"

        idx += count
    nat[nat == "TS"] = ""

    # Extratropical transition
    # Look after the last TC point for ET
    if et:
        idx = np.where(nat == "TC")[0]
        if len(idx) > 0:
            idx = idx[-1] + 1
            new_idx = fill_next_nature(
                nat, ~symmetric & warm_core & ~weak, "ET", idx, min_count
            )
            if new_idx >= idx + min_count:
                idx = new_idx
                new_idx = fill_next_nature(
                    nat, ~symmetric & ~warm_core & ~weak, "BC", idx, min_count
                )
                if new_idx >= idx + min_count:
                    idx = new_idx
                    fill_next_nature(nat, warm_core & ~weak, "WS", idx, min_count)
                else:
                    # If nothing was labelled as baroclinic following ET, remove ET
                    nat[np.isin(nat, ["ET", "BC"])] = ""
            else:
                # ET lasted less than min_count remove ET
                nat[nat == "ET"] = ""

    # Label remaining unlabelled sections
    # Asymmetric = baroclinic
    nat[(nat == "") & ~symmetric] = "BC"
    # Upper-level cold core = Trough
    nat[(nat == "") & trough] = "Tr"
    # Low-level cold core = Mid level vortex
    nat[(nat == "") & ~warm_core] = "MV"
    # Warm core symmetric not TC (decaying)
    nat[(nat == "")] = "Ot"

    if smooth:
        smooth_excursions(nat, min_count)

    return nat

