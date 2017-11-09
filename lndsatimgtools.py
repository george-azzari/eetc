import ee


def _rename_band(val, suffix):
    return ee.String(val).cat(ee.String("_")).cat(ee.String(suffix))


def rename_bands(img, suffix):
    bandnames = img.bandNames()
    newnames = bandnames.map(lambda x: _rename_band(x, suffix))
    return img.select(bandnames, newnames)


def add_latlon(img):
    ll = img.select('BLUE').multiply(0).add(ee.Image.pixelLonLat())
    return img.addBands(ll.select(['longitude', 'latitude'], ['LON', 'LAT']))


def add_gcviband(img):
    b2 = img.select('GREEN').toFloat()
    b4 = img.select('NIR').toFloat()
    gcvi = (b4.divide(b2)).subtract(1.0)
    return img.addBands(gcvi.select([0], ['GCVI']))


def add_logcviband(img):
    b2 = img.select('GREEN').toFloat()
    b4 = img.select('NIR').toFloat()
    gcvi = (b4.divide(b2)).subtract(1.0).log()
    return img.addBands(gcvi.select([0], ['LOGGCVI']))


def add_ndviband(img):
    ndvi = img.normalizedDifference(['NIR', 'RED'])
    return img.addBands(ndvi.select([0], ['NDVI']))


def count_region(img, region, scale=30):
    reduced = img.select([0], ['count']).reduceRegion(
        reducer=ee.Reducer.count(),
        geometry=region,
        scale=scale,
        maxPixels=1e13)
    return ee.Number(reduced.get('count'))
