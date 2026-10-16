# Glacier snow cover fraction from Sentinel-1, Sentinel-2 and MODIS

Python and Jupyter implementation of the multi-sensor glacier snow-covered area fraction (SCAF) workflow described in:

> Kim, D., Mattea, E., Callegari, M., Saks, T., Kenzhebayev, R., Azisov, E., Ullmann, T., Hoelzle, M., and Barandun, M.: **Sub-seasonal snowline dynamics of glaciers in Central Asia from multi-sensor satellite observations, 2000–2023**, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-3978, 2025.

The workflow derives glacier SCAF from Sentinel-1 radar and Sentinel-2 optical observations, then uses these estimates to calibrate an empirical relationship between MODIS mean near-infrared (NIR) reflectance and SCAF. Applying that relationship to the longer MODIS record supports high-frequency monitoring of snow depletion on glaciers.

Satellite processing runs in **Google Earth Engine (GEE)**. Time series are returned as date-indexed `pandas.Series`; regression fitting and plotting run locally.

## Repository contents

| File | Purpose |
| --- | --- |
| [`glacier_scaf.py`](glacier_scaf.py) | Earth Engine extraction, Sentinel classification and snowline correction, and MODIS–Sentinel regression fitting. |
| [`timeseries_extract.ipynb`](timeseries_extract.ipynb) | Extract and plot MODIS mean NIR and Sentinel-1/-2 SCAF for one glacier and a specific period. |
| [`MODIS_SCAF_Sentinel_regression.ipynb`](MODIS_SCAF_Sentinel_regression.ipynb) | Extract MODIS mean NIR and Sentinel-1 and Sentinel-2 SCF for 2016–2023 and fit the MODIS NIR–SCAF relationship. |
| [`MODIS_SCAF.ipynb`](MODIS_SCAF.ipynb) | Convert MODIS mean NIR to SCAF using specified regression coefficients. |


## Setup

You need Python, a Jupyter environment, an Earth Engine-enabled account and Google Cloud project, and access to a glacier outline in Earth Engine. Keep `glacier_scaf.py` and the notebooks in the same directory, and run Jupyter from that directory.

Create an environment and install the packages used by the code and notebooks:

```bash
python -m venv .venv
source .venv/bin/activate
# Windows: .venv\Scripts\activate
python -m pip install earthengine-api geemap pandas numpy scipy matplotlib jupyterlab
jupyter lab
```

Authenticate and initialize Earth Engine **before importing `glacier_scaf`**:

```python
import ee

ee.Authenticate()  # Run when credentials need to be established or renewed.
ee.Initialize(project="YOUR_EARTH_ENGINE_PROJECT")

from glacier_scaf import MODIS_meanNIR, SCAF_S1, SCAF_S2, fit_scaf_regression
```

Replace `YOUR_EARTH_ENGINE_PROJECT` in each notebook with your project ID. The project must be configured for Earth Engine and your account must have permission to use it. See the [Earth Engine Python installation and authentication guide](https://developers.google.com/earth-engine/guides/python_install).

## Function reference

| Function | Return value |
| --- | --- |
| `MODIS_meanNIR(glacier_outline, date_start, date_end)` | Sorted, date-indexed `pandas.Series` of scaled mean NIR reflectance, with rejected values represented as `NaN`. |
| `SCAF_S1(start_date, end_date, glacier_outline, orbit)` | Per-acquisition, date-indexed SCAF percentages, named `scaf_S1`. |
| `SCAF_S2(start_date, end_date, glacier_outline, th_min, th_max, th_fixed, tile)` | Per-acquisition, date-indexed SCAF percentages, named `scaf_S2`. |
| `fit_scaf_regression(modis_nir, s1_scaf, s2_scaf)` | Tuple of Python floats `(a, b)`. Fitting itself makes no Earth Engine requests. |
