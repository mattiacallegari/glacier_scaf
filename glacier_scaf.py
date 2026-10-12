"""Google Earth Engine glacier NIR and snow cover time series.

The caller must authenticate and initialize Earth Engine before calling these
functions. Dates are inclusive; time series values are per acquisition (0–100
percent for SCAF). Maps use 1 = snow, 0 = ice, masked = no observation.
"""

from datetime import date, datetime, timedelta

import ee
import pandas as pd


DEM = "NASA/NASADEM_HGT/001"
SLA_BIN_M = 25


def _dates(start_date, end_date):
    start = pd.Timestamp(start_date).date()
    end = pd.Timestamp(end_date).date()
    if end < start:
        raise ValueError("end_date must be on or after start_date")
    return start.isoformat(), (end + timedelta(days=1)).isoformat()


def _geometry(outline):
    if isinstance(outline, ee.geometry.Geometry):
        return outline
    if isinstance(outline, (ee.feature.Feature, ee.featurecollection.FeatureCollection)):
        return outline.geometry()
    raise TypeError("glacier_outline must be an ee.Geometry, ee.Feature or ee.FeatureCollection")


def _number_or(value, default):
    return ee.Number(ee.Algorithms.If(value, value, default))


def _fraction(image, geometry, scale):
    """Percent snow among observed snow and ice pixels, or null when empty."""
    sums = image.select("classification").reduceRegion(
        reducer=ee.Reducer.sum().combine(ee.Reducer.count(), sharedInputs=True),
        geometry=geometry, scale=scale, maxPixels=1_000_000_000,
    )
    count = _number_or(sums.get("classification_count"), 0)
    snow = _number_or(sums.get("classification_sum"), 0)
    return ee.Algorithms.If(count.gt(0), snow.divide(count).multiply(100), None)


def _series(collection, property_name, name):
    """Fetch the small results table once; retain multiple acquisitions per day."""
    def row(image):
        image = ee.Image(image)
        return ee.Feature(None, {
            "date": ee.Date(image.get("system:time_start")).format("YYYY-MM-dd'T'HH:mm:ss"),
            "value": image.get(property_name),
        })

    rows = ee.FeatureCollection(collection.toList(collection.size()).map(row)).getInfo()["features"]
    data = [(f["properties"]["date"], f["properties"]["value"])
            for f in rows if f["properties"].get("value") is not None]
    if not data:
        return pd.Series([], index=pd.DatetimeIndex([], name="date"),
                         dtype="float64", name=name)
    series = pd.Series([float(v) for _, v in data],
                       index=pd.DatetimeIndex([d for d, _ in data], name="date"),
                       name=name, dtype="float64")
    return series.sort_index()


def MODIS_meanNIR(start_date, end_date, glacier_outline):
    """Daily MOD09GQ band 2 mean NIR (native scaled integers, approximately 0–10000).

    MOD09GA state_1km cloud-state bits mask cloudy (state 1) and missing
    (state 3) pixels; only dates with at least 75% valid 250 m glacier area
    are kept. A 250 m pixel is included if its center lies in the outline.
    """
    start, stop = _dates(start_date, end_date)
    geom = _geometry(glacier_outline)
    nir = ee.ImageCollection("MODIS/061/MOD09GQ").filterDate(start, stop).filterBounds(geom)
    qa = ee.ImageCollection("MODIS/061/MOD09GA").filterDate(start, stop).filterBounds(geom)
    paired = ee.Join.saveFirst("qa_image").apply(
        nir, qa, ee.Filter.equals(leftField="system:time_start", rightField="system:time_start")
    )
    paired = ee.ImageCollection(paired).filter(ee.Filter.notNull(["qa_image"]))

    def measure(image):
        image = ee.Image(image)
        state = ee.Image(image.get("qa_image")).select("state_1km").bitwiseAnd(3)
        valid = state.neq(1).And(state.neq(3))
        band = image.select("sur_refl_b02").updateMask(valid).updateMask(
            image.select("sur_refl_b02").gte(0))
        stats = band.reduceRegion(
            reducer=ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True),
            geometry=geom, scale=250, maxPixels=1_000_000_000,
        )
        total = ee.Image.constant(1).rename("area").reduceRegion(
            reducer=ee.Reducer.count(), geometry=geom, scale=250,
            crs=image.select("sur_refl_b02").projection(), maxPixels=1_000_000_000,
        ).get("area")
        count = _number_or(stats.get("sur_refl_b02_count"), 0)
        enough = count.gte(ee.Number(total).multiply(0.75))
        return image.set("meanNIR", ee.Algorithms.If(enough,
                         stats.get("sur_refl_b02_mean"), None))

    return _series(paired.map(measure), "meanNIR", "MODIS_meanNIR")


def _otsu(histogram, fallback):
    """Server-side Otsu threshold from an EE histogram dictionary."""
    hist = ee.Dictionary(histogram)
    counts = ee.Array(hist.get("histogram"))
    means = ee.Array(hist.get("bucketMeans"))
    size = ee.Number(means.length().get([0]))
    total = ee.Number(counts.reduce(ee.Reducer.sum(), [0]).get([0]))
    weighted = ee.Number(means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]))

    def score(i):
        i = ee.Number(i)
        left = counts.slice(0, 0, i)
        n = ee.Number(left.reduce(ee.Reducer.sum(), [0]).get([0]))
        right_n = total.subtract(n)
        left_sum = ee.Number(means.slice(0, 0, i).multiply(left)
                             .reduce(ee.Reducer.sum(), [0]).get([0]))
        left_mean = left_sum.divide(n.max(1))
        right_mean = weighted.subtract(left_sum).divide(right_n.max(1))
        return ee.Number(ee.Algorithms.If(
            n.gt(0).And(right_n.gt(0)),
            n.multiply(left_mean.subtract(weighted.divide(total)).pow(2))
             .add(right_n.multiply(right_mean.subtract(weighted.divide(total)).pow(2))),
            -1))

    indices = ee.List.sequence(1, size.subtract(1))
    best = ee.Number(indices.sort(indices.map(score)).get(-1))
    return ee.Number(ee.Algorithms.If(size.gt(1), means.get([best]), fallback))


def _topographic_b8(image, geom, dem):
    """SCS+C correction of S2 B8, falling back to B8 if fit is undefined."""
    radians = ee.Number(3.141592653589793).divide(180)
    zenith = ee.Number(image.get("MEAN_SOLAR_ZENITH_ANGLE")).multiply(radians)
    azimuth = ee.Number(image.get("MEAN_SOLAR_AZIMUTH_ANGLE")).multiply(radians)
    slope = ee.Terrain.slope(dem).multiply(radians)
    aspect = ee.Terrain.aspect(dem).multiply(radians)
    cos_z = zenith.cos()
    cos_s = slope.cos()
    illumination = cos_s.multiply(cos_z).add(
        slope.sin().multiply(zenith.sin()).multiply(azimuth.subtract(aspect).cos()))
    b8 = image.select("B8")
    fit = illumination.rename("IC").addBands(b8).updateMask(
        slope.gte(ee.Number(5).multiply(radians)).And(illumination.gte(0))
    ).reduceRegion(reducer=ee.Reducer.linearFit(), geometry=geom.buffer(1000),
                   scale=10, maxPixels=1_000_000_000, bestEffort=True)
    slope_fit = _number_or(fit.get("scale"), 0)
    c = _number_or(fit.get("offset"), 0).divide(
        ee.Number(ee.Algorithms.If(slope_fit.abs().gt(1e-8), slope_fit, 1)))
    corrected = b8.multiply(cos_s.multiply(cos_z).add(c)).divide(illumination.add(c))
    use_fit = slope_fit.abs().gt(1e-8)
    return ee.Image(ee.Algorithms.If(use_fit, corrected, b8)).rename("B8") \
        .unmask(b8).copyProperties(image, ["system:time_start", "system:index"])


def SLA_correction(SC_map):
    """Correct an S1 or S2 EE snow map carrying glacier_geometry and sensor.

    Estimate the first altitude bin with >=50% snow in three consecutive
    25 m bins. Beyond 20% of glacier elevation range above/below this line,
    flip implausible labels; S2 follows the source notebook's upper-only rule
    and applies it only when preliminary SCAF >=25%. If no line is found,
    leave the image unchanged and set its 'sla' property to -9999.
    """
    image = ee.Image(SC_map)
    geom = ee.Geometry(image.get("glacier_geometry"))
    dem = ee.Image(DEM).select("elevation")
    bounds = dem.reduceRegion(reducer=ee.Reducer.minMax(), geometry=geom,
                              scale=30, maxPixels=1_000_000_000)
    low = ee.Number(bounds.get("elevation_min"))
    high = ee.Number(bounds.get("elevation_max"))
    start = low.divide(SLA_BIN_M).floor().multiply(SLA_BIN_M)
    count = high.subtract(start).divide(SLA_BIN_M).ceil().max(0)
    levels = ee.List.sequence(0, count)
    classification = image.select("classification")

    def share(i):
        base = start.add(ee.Number(i).multiply(SLA_BIN_M))
        mask = dem.gte(base).And(dem.lt(base.add(SLA_BIN_M)))
        stats = classification.updateMask(mask).reduceRegion(
            reducer=ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True),
            geometry=geom, scale=ee.Number(image.get("analysis_scale")),
            maxPixels=1_000_000_000)
        n = _number_or(stats.get("classification_count"), 0)
        return ee.Number(ee.Algorithms.If(n.gt(0), stats.get("classification_mean"), -1))

    shares = levels.map(share)
    candidates = ee.List.sequence(0, count.subtract(2).max(0)).map(
        lambda i: ee.Algorithms.If(
            ee.Number(shares.get(ee.Number(i))).gte(0.5)
            .And(ee.Number(shares.get(ee.Number(i).add(1))).gte(0.5))
            .And(ee.Number(shares.get(ee.Number(i).add(2))).gte(0.5)), i, -1)
    ).removeAll([-1])
    found = candidates.size().gt(0).And(count.gte(2))
    sla = ee.Number(ee.Algorithms.If(found,
                    start.add(ee.Number(candidates.get(0)).multiply(SLA_BIN_M)), -9999))
    margin = high.subtract(low).multiply(0.2)
    upper = classification.eq(0).And(dem.gte(sla.add(margin))).And(found)
    lower = classification.eq(1).And(dem.lt(sla.subtract(margin))).And(found)
    is_s2 = ee.String(image.get("sensor")).compareTo("S2").eq(0)
    s2_allowed = ee.Number(_number_or(image.get("scf_glacier"), 0)).gte(25)
    upper = upper.And(ee.Image.constant(ee.Number(ee.Algorithms.If(is_s2, s2_allowed, 1))))
    lower = lower.And(ee.Image.constant(ee.Number(ee.Algorithms.If(is_s2, 0, 1))))
    result = classification.where(upper, 1).where(lower, 0).clip(geom)
    return result.copyProperties(image, image.propertyNames()).set({
        "sla": sla,
        "scafAfter": _fraction(result, geom, ee.Number(image.get("analysis_scale"))),
    })


def SCAF_S2(start_date, end_date, glacier_outline, th_min, th_max, th_fixed,
            SLA_filetr=True):
    """Cloud screened S2 L1C B8 Otsu snow fraction (percent per acquisition).

    Cloud probability >40 is cloud; discard scenes with >10% cloudy pixels
    on the glacier. Threshold falls back to th_fixed outside (th_min, th_max).
    """
    if not th_min < th_max or not th_min < th_fixed < th_max:
        raise ValueError("Expected th_min < th_fixed < th_max (S2 B8 digital numbers)")
    start, stop = _dates(start_date, end_date)
    geom = _geometry(glacier_outline)
    s2 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED").filterDate(start, stop).filterBounds(geom)
    cloud = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY").filterDate(start, stop).filterBounds(geom)
    joined = ee.ImageCollection(ee.Join.saveFirst("cloud_image").apply(
        s2, cloud, ee.Filter.equals(leftField="system:index", rightField="system:index")
    )).filter(ee.Filter.notNull(["cloud_image"]))
    dem = ee.Image(DEM).select("elevation")

    def classify(image):
        image = ee.Image(image)
        clouds = ee.Image(image.get("cloud_image")).select("probability").gt(40)
        cloud_cover = clouds.reduceRegion(reducer=ee.Reducer.mean(), geometry=geom,
                                          scale=10, maxPixels=1_000_000_000).get("probability")
        nir = _topographic_b8(image, geom, dem).updateMask(clouds.Not()).clip(geom)
        histogram = nir.reduceRegion(reducer=ee.Reducer.histogram(50, 2), geometry=geom,
                                     scale=10, maxPixels=1_000_000_000).get("B8")
        threshold = ee.Number(ee.Algorithms.If(histogram, _otsu(histogram, th_fixed), th_fixed))
        threshold = ee.Number(ee.Algorithms.If(threshold.gt(th_min).And(threshold.lt(th_max)),
                                               threshold, th_fixed))
        snow = nir.gte(threshold).rename("classification")
        return snow.set({"cloud_percent_glacier": _number_or(cloud_cover, 100).multiply(100),
                         "otsu_threshold": threshold, "scf_glacier": _fraction(snow, geom, 10),
                         "sensor": "S2", "analysis_scale": 10,
                         "glacier_geometry": geom}) \
            .copyProperties(image, ["system:time_start", "system:index"])

    mapped = joined.map(classify).filter(ee.Filter.lte("cloud_percent_glacier", 10))
    if SLA_filetr:
        mapped = mapped.map(SLA_correction)
    return _series(mapped, "scafAfter" if SLA_filetr else "scf_glacier", "SCAF_S2")


def SCAF_S1(start_date, end_date, glacier_outline, SLA_filetr=True):
    """S1 IW ascending VV/VH wet snow fraction (percent per acquisition).

    Use the same relative orbit's images within 30 days of January 5 in
    each acquisition year as a winter reference. Mean VV/VH differences below
    -6 dB are wet snow; classify at 50 m. Acquisitions without a winter
    reference for their relative orbit are omitted.
    """
    start, stop = _dates(start_date, end_date)
    geom = _geometry(glacier_outline)
    base = ee.ImageCollection("COPERNICUS/S1_GRD").filterBounds(geom) \
        .filter(ee.Filter.eq("instrumentMode", "IW")) \
        .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV")) \
        .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH")) \
        .filter(ee.Filter.eq("orbitProperties_pass", "ASCENDING"))
    scenes = base.filterDate(start, stop)

    def classify(image):
        image = ee.Image(image)
        ref_day = ee.Date.fromYMD(ee.Date(image.get("system:time_start")).get("year"), 1, 5)
        reference = base.filterDate(ref_day.advance(-30, "day"),
                                    ref_day.advance(30, "day")) \
            .filter(ee.Filter.eq("relativeOrbitNumber_start",
                                  image.get("relativeOrbitNumber_start")))
        # Conditional EE image prevents empty-reference reductions from failing.
        winter = ee.Image(ee.Algorithms.If(reference.size().gt(0), reference.select(
            ["VV", "VH"]).mean(), image.select(["VV", "VH"])))
        ratio = image.select("VV").subtract(winter.select("VV")) \
            .add(image.select("VH").subtract(winter.select("VH"))).divide(2)
        coarse = ratio.reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024) \
            .reproject(crs=image.select("VV").projection(), scale=50)
        snow = coarse.lt(-6).rename("classification").clip(geom)
        return snow.set({"reference_count": reference.size(),
                         "scf_glacier": _fraction(snow, geom, 50),
                         "sensor": "S1", "analysis_scale": 50,
                         "glacier_geometry": geom}) \
            .copyProperties(image, ["system:time_start", "system:index"])

    mapped = scenes.map(classify).filter(ee.Filter.gt("reference_count", 0))
    if SLA_filetr:
        mapped = mapped.map(SLA_correction)
    return _series(mapped, "scafAfter" if SLA_filetr else "scf_glacier", "SCAF_S1")
