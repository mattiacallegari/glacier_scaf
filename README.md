# Glacier NIR and snow cover time series

This repository packages the methods in the supplied MODIS, Sentinel-1, and Sentinel-2 notebooks as reusable Python functions. Calculations run on Google Earth Engine (EE); the functions return local `pandas.Series` indexed by UTC acquisition timestamps. The example notebook plots the 2019 melt season for Vernagtferner.

## Setup

Python 3.10+ and an Earth Engine account with access to a Cloud project are required. Install dependencies:

```bash
python -m pip install earthengine-api pandas matplotlib jupyter
```

Authenticate once, then initialize the EE client in each new Python session:

```python
import ee
# Run ee.Authenticate() once if you have not authenticated on this machine.
ee.Initialize(project="YOUR_EARTH_ENGINE_PROJECT")
```

Run `jupyter notebook example_2019.ipynb` from this directory and replace the project placeholder. Internet access and EE permissions are needed to execute the notebook; calculations may take several minutes and very long periods may exceed EE computation or response quotas.

## API

```python
from glacier_timeseries import MODIS_meanNIR, SCAF_S2, SCAF_S1, SLA_correction

outline = ee.FeatureCollection("GLIMS/current").filter(
    ee.Filter.eq("glac_id", "G010823E46875N")
).sort("src_date", False).first().geometry()

nir = MODIS_meanNIR("2019-06-01", "2019-09-30", outline)
s2 = SCAF_S2("2019-06-01", "2019-09-30", outline,
             th_min=3500, th_max=5000, th_fixed=4500)
s1 = SCAF_S1("2019-06-01", "2019-09-30", outline)
```

| Function | Result | Method |
| --- | --- | --- |
| `MODIS_meanNIR(start_date, end_date, glacier_outline)` | Series of MOD09GQ band 2 native reflectance integers (roughly 0–10000) | MOD09GA 1 km QA cloud state excludes cloud and missing data; requires ≥75% valid 250 m glacier pixels. |
| `SCAF_S2(start_date, end_date, glacier_outline, th_min, th_max, th_fixed, SLA_filetr=True)` | Series of snow fraction in percent (0–100) | Sentinel-2 L1C harmonized B8, cloud probability ≤40 and ≤10% clouds over glacier, SCS+C terrain correction, per-scene Otsu threshold with fixed fallback. Thresholds use native B8 digital numbers. |
| `SCAF_S1(start_date, end_date, glacier_outline, SLA_filetr=True)` | Series of wet snow fraction in percent (0–100) | Sentinel-1 ascending IW VV/VH, each orbit compared with its January winter reference; average difference <−6 dB at 50 m indicates wet snow. |
| `SLA_correction(SC_map)` | `ee.Image` with corrected `classification` band and `sla`, `scafAfter` properties | Shared elevation bin algorithm. See below for input requirements. |

`glacier_outline` accepts an `ee.Geometry`, `ee.Feature`, or `ee.FeatureCollection`. Dates are **inclusive**. The series include acquisition times, so multiple same-day acquisitions remain distinct. Missing measurements and S1 scenes with no winter reference are omitted; the output is an empty float Series when no scenes qualify. `SLA_filetr` is spelled as requested in the API; set it to `False` to omit SLA correction.

For a direct call to `SLA_correction`, supply an `ee.Image` with a `classification` band (snow=1, ice=0, missing pixels masked) and properties `glacier_geometry` (`ee.Geometry`), `sensor` (`"S1"` or `"S2"`), and `analysis_scale` (metres). For S2 also supply `scf_glacier` as a percentage. The result is lazy: call `getInfo()` or export it to evaluate it. An absent qualifying snowline leaves the classification unchanged and sets `sla=-9999`.

SLA follows three consecutive 25 m altitude bins containing at least 50% snow. Above a margin of 20% of the glacier elevation range, S1 and S2 can reclassify ice as snow. Below that margin S1 can also reclassify snow as ice. The S2 upper correction runs only when uncorrected SCAF is at least 25%, following its source notebook.

## Notes on the source notebooks

The notebooks include hardcoded tiles, glacier IDs, dates, export paths, and exploratory cells. This module takes a glacier geometry and dates instead and performs no exports. It preserves the core datasets and classification steps. MODIS cloud screening uses a valid-pixel fraction in place of the notebook's pixel-grid selection and one-cloud-pixel rule; S1 uses January 5 in each scene's year rather than the notebook's fixed January 2019 reference. The S2 scene filter accepts any intersecting MGRS tile, making the API usable beyond the source example. Scenes may repeat across tiles or orbits.

The original notebooks are provided separately as reference material. This example has not been run against a live Earth Engine project here; cloud availability, GLIMS collection access, and computational limits depend on your EE account.
