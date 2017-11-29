import ee

# def addTerrain(image):
#
#      terrain = ned.clip(geometry).reproject(image.projection());
#
#      provterrain = terrain.updateMask(image.select(0).mask()).select([0], ['elev']);
#      provaspect = ee.Terrain.aspect(provterrain).select([0], ['aspect']);
#      provslope = ee.Terrain.slope(provterrain).select([0], ['slope']);
#
#     return composite
#         .addBands(provterrain)
#         .addBands(provaspect)
#         .addBands(provslope)


# Add two bands represeting lon/lat of each pixels
def add_latlon(image):

    ll = image.select(0).multiply(0).add(ee.Image.pixelLonLat())

    return image.addBands(ll.select(['longitude', 'latitude'], ['LON', 'LAT']))


# TODO: this is incomplete, needs review and to be extended to custom classes
def stratify_cdl(year, npoints, region, allcrops, random_seed=1234):

    cdlcoll = ee.ImageCollection("USDA/NASS/CDL")

    cdl = ee.Image(cdlcoll.filter(ee.Filter.calendarRange(year, year, 'YEAR')).first())
    cdl = cdl.select(['cropland'], ['CDL'])

    cropmask = cdl.gte(196).add(cdl.lte(60))
    cropmask = cropmask.add(cdl.gte(66).multiply(cdl.lte(77))).select([0], ['CROPMASK'])

    tosample = add_latlon(cdl.updateMask(cropmask).addBands(cropmask))

    samples = tosample.stratifiedSample(
        numPoints=npoints,
        classBand='CROPMASK',
        region=region.geometry(),
        scale=30,
        seed=random_seed,
        classValues=[0, 1],
        classPoints=[0, npoints],
        tileScale=16)
