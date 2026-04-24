# Extended IBTrACS

Author: Stella Bourdin (Uni. Oxford), with help from Leo Saffin and Kevin Hodges (Uni. Reading), and advice from many members of the [Huracán project](https://research.reading.ac.uk/huracan/about-huracan/project-overview/).

This code and data are distributed under the GPL-3.0 license. This does not cover the IBTrACs data, whose license is not indicated on the website. 

A paper is being written to be submitted to ESSD describing the construction process of the present database.

## About Git LFS
Due to the size of the files in this repository, we use git Large File Storage (LFS) to store the input and output files. 
To retrieve these files from this repository, you will need to install [git LFS](https://git-lfs.com/).
Then, once you cloned the repository, run `git lfs pull` to download all the files (it will take time). 

## Repository description
* `ibtracs/` contains csv files of IBTrACS data for each basin;
* `input/` contains pickle files of track data for each source;
* `extended-ibtracs/` contains NetCDF files of the extended-ibtracs data for each basin;
* `eib.py` is a small python module to quickly import the packages used in constructing and analyse the data;
* `*.ipynb` are the Jupyter Notebooks used in constructing the database. 

## Data description

Coordinates:
* `source` is the origin of the track data, either IBTrACS, or a (re)-analysis.
* `record` is an index for each TC point across all sources.

Data variables:
* `track_id` is the unique TC identifier from IBTrACS.
* `time` is the time of the corresponding record.
* `lon` and `lat` describe the position of the point in `track_id` at `time` for each source
* `pres` and `wind10` describe the intensity (in terms of minimum SLP and maximum surface wind speed, respectively) of the point in `track_id` at `time` for each source. NB : this variable is empty for the IBTrACS source, where there are several wind attributes, and it is up to the user to decide which one they may use.
* `translation_speed` and `azimuth` have been computed along the tracks; they are in m/s and degrees, respectively.
* `vtu`, `vtl` and `b` are the Cyclone Phase Space (CPS, Hart 2003) parameters. They are currently only provided for some of the TRACK sources.
* `vtu_roll`, `vtl_roll` and `b_roll` are the 24h rolling means of the CPS parameters.
* `WCSI` means Warm Core Symmetric Intensifying. It is True when the cyclone has a Symmetric Warm Core as per the CPS parameters (|B|<15 & VTL>0), and when it is intensifying at the same time (intensification being defined from the 850hPa vorticity)
* `CCA` means Cold Core Assymetric, and is True when the cyclone has an Assymetric Cold Core as per the CPS parameters (|B|>15 & VTL<0).
* `short_label` is the SyCLoPS label for Low Pressure System category. It is only provided for SyCLoPS sources.
* `is_tc` is a composite flag for identifying the TC stage. For IBTrACS, it is True when `nature=='TS'`, for SyCLoPS, it is True when `short_label=='TC'`, and for TRACK it is True when WCSI is True.
* `is_etc` is a composite flag for identifying the TC stage. For IBTrACS, it is True when `nature=='ET'`, for SyCLoPS, it is True when `short_label=='EX'`, and for TRACK it is True when CCA is True.
* `ET` flags the onset and offset of Extra-Tropical Transition. `ET=+1` corresponds to the onset of transition, and corresponds to the last timestep when `is_tc` is True for a given track_id and source. `ET=-1` corresponds to the offset of transition, and corresponds to the first timestep when `is_etc` is True for a given track_id and source. 
* All other attributes along the dimension `record` only come from IBTrACS, so please refer to the [IBTrACS Column Documentation](https://www.ncei.noaa.gov/sites/default/files/2021-07/IBTrACS_v04_column_documentation.pdf).
