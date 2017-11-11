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
