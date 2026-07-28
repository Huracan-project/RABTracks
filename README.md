# RABTracks

Author: Stella Bourdin (Formerly Uni. Oxford, now LMD-IPSL, Ecole Polytechnique), with help from Leo Saffin and Kevin Hodges (Uni. Reading), and Yushan Han (UCDavis), and advice from many members of the [Huracán project](https://research.reading.ac.uk/huracan/about-huracan/project-overview/).

This code and data are distributed under the GPL-3.0 license. This does not cover the IBTrACs data, whose license is not indicated on the website. 

A paper is being written to be submitted to ESSD describing the construction process of the present database. Pre-print is available on [HAL](https://hal.science/hal-05687753v1).

## Obtaining the RABTracks dataset

The RABTracks dataset is publicly available through [this link](https://gws-access.jasmin.ac.uk/public/huracan/RABTracks/). Download the `rabtracks.zip` archive.

The archive contains one file per basin. Each file has been compressed using the LZMA algorithm. They can be decompressed using the `unxz` or `xz --decompress` command. 

### Data description

**Coordinates**:
* `source` is the origin of the track data, either IBTrACS, or a (re)-analysis.
* `record` is an index for each TC point across all sources.

**Data variables**:
* `track_id` is the unique TC identifier from IBTrACS.
* `time` is the time of the corresponding record.
* `lon` and `lat` describe the position of the point in `track_id` at `time` for each source. `lon_pres` & `lat_pres` and `lon_wind` & `lat_wind` describe respectively the position of the pressure minimum and wind maximum associated with the TC center. 
* `pres` and `wind10` describe the intensity (in terms of minimum SLP and maximum surface wind speed, respectively) of the point in `track_id` at `time` for each source. NB : this variable is empty for the IBTrACS source, where there are several wind attributes, and it is up to the user to decide which one they may use.
* `vo850` describes the intensity in terms of relative vorticity at 850hPa.
* `translation_speed` and `azimuth` have been computed along the tracks; they are in m/s and degrees, respectively. The `*_roll` counterparts have been computed applying a 24h rolling mean to the lon and lat data.
* `radius_max_wind` is the radius of maximum wind, computed as the distance between the pressure minimum and the wind maximum, where these position were available (for TRACK data). The `radius_max_wind_roll` counterpart as been smoothed with a 24h rolling mean.
* `lpsarea` is the Low Pressure System area as computed by SyCLoPS.
* `vtu`, `vtl` and `b` are the Cyclone Phase Space (CPS, Hart 2003) parameters. They are currently only provided for some of the TRACK sources.
* `vtu_roll`, `vtl_roll` and `b_roll` are the 24h rolling means of the CPS parameters.
* `WCSI` means Warm Core Symmetric Intensifying. It is True when the cyclone has a Symmetric Warm Core as per the CPS parameters (|B|<15 & VTL>0), and when it is intensifying at the same time (intensification being defined from the 850hPa vorticity)
* `CCA` means Cold Core Assymetric, and is True when the cyclone has an Assymetric Cold Core as per the CPS parameters (|B|>15 & VTL<0).
* `label` is the SyCLoPS (adjusted) label for Low Pressure System category. It is only provided for SyCLoPS sources.
* `is_tc` is a composite flag for identifying the TC stage. For IBTrACS, it is True when `nature=='TS'`, for SyCLoPS, it is True when `short_label=='TC'`, and for TRACK it is True when WCSI is True (see ESSD paper for more details).
* `is_etc` is a composite flag for identifying the TC stage. For IBTrACS, it is True when `nature=='ET'`, for SyCLoPS, it is True when `short_label=='EX'`, and for TRACK it is True when CCA is True (see ESSD paper for more details).
* `ET` flags the onset and offset of Extra-Tropical Transition. `ET=+1` corresponds to the onset of transition, and corresponds to the last timestep when `is_tc` is True for a given track_id and source. `ET=-1` corresponds to the offset of transition, and corresponds to the first timestep when `is_etc` is True for a given track_id and source  (see ESSD paper for more details). 
* All other attributes along the dimension `record` only come from IBTrACS, so please refer to the [IBTrACS Column Documentation](https://www.ncei.noaa.gov/sites/default/files/2021-07/IBTrACS_v04_column_documentation.pdf).

### Examples of how to read/use
You can check the notebooks in `paper_figs` to see the code used to generate the figures for the ESSD paper. 

## Building RABTracks yourself
If you want to run the code yourself, **(1) clone this repository**, and then **(2) download the input data**:
it is available through the [same link](https://gws-access.jasmin.ac.uk/public/huracan/RABTracks/), but you will need to download the `input.zip` archive.

Decompress the zip file into your RABTracks folder, and decompress all the compressed `.xz` files using the `unxz` or `xz --decompress` command. If necessary, move the `ibtracs` and `RA_tracks` folders to the first level, and create the `RABTracks` folder to store the output. Your folder structure should be similar to the one below. From there, you can run `2_build_RABTracks.ipynb` (and modify it if you like) to build the dataset. Note that you might need to install a python environment with the necessary packages.

```
.
├── ibtracs -> Contains the input files from IBTrACS
│   ├── ibtracs_EP.pkl
│   ├── ibtracs_EP.pkl.xz
│   ├── ibtracs_NA.pkl
│   ├── ibtracs_NA.pkl.xz
│   ├── ...
├── RA_tracks -> Contains the input files from the different reanalyses/trackers
│   ├── SyCLoPS-ERA5.pkl
│   ├── SyCLoPS-ERA5.pkl.xz
│   ├── ...
│   ├── TRACK_CPS-ERA5.pkl
│   ├── TRACK_CPS-ERA5.pkl.xz
│   ├── ...
│   ├── TRACK-ECMWF-OP-AN.pkl
│   ├── TRACK-ECMWF-OP-AN.pkl.xz
│   ├── ...
│   ├── TRACK_WCSI-ERA5.pkl
│   ├── TRACK_WCSI-ERA5.pkl.xz
│   ├── ...
├── RABTracks -> Will store the output
│   ├── Create this folder if necessary!
├── paper_figs -> Script to create the figures (examples of how to read/use the data)
│   ├── ...
├── LICENSE
├── rabtracks.py
├── 1_input_data.ipynb -> Was used to create the input from raw data (not published)
├── 2_build_RABTracks.ipynb -> Use to create RABTracks from the input
├── 2_update_public.txt
└── README.md
```

## References
SyCLoPS:
- Paper [Han & Ullrich, JGR: The System for Classification of Low-Pressure Systems (SyCLoPS): An All-In-One Objective Framework for Large-Scale Data Sets](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024JD041287)
- [GitHub](https://github.com/yepkids/SyCLoPS)
