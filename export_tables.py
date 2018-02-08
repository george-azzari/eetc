"""
Author: George Azzari (gazzari@stanford.edu)
Center on Food Security and the Environment
Department of Earth System Science
Stanford University
"""

import ee


def export_features(features, fname, export_to='drive'):

    if export_to == 'drive':
        task = ee.batch.Export.table.toDrive(features, fname, '')

    else:
        task = ee.batch.Export.table.toCloudStorage(features,
                                                    description=fname,
                                                    bucket='us-cdl-samples',
                                                    fileNamePrefix=fname,
                                                    fileFormat=None)
    task.start()

    return task


def get_point(image, point, scale, depth):

    s = image.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=point.geometry(),
        scale=scale,
        bestEffort=False,
        maxPixels=depth,
        tileScale=16)

    return ee.Feature(ee.Feature(None, s).copyProperties(image).copyProperties(point)).set({

        'PTLON': ee.Number(ee.List(point.geometry().coordinates()).get(0)),
        'PTLAT': ee.Number(ee.List(point.geometry().coordinates()).get(1)),
        'OBSDATE': ee.Date(image.get('system:time_start')).format(),
        'MSTIME': image.get('system:time_start')

    })


def sampleGridCore(image, grid, scale):

    depth = image.bandNames().length()
    samples = grid.map(lambda point: get_point(image, point, scale, depth))

    return samples


def sampleGridTS(imagecoll, grid, scale, control, doexport, fname):

    samples = imagecoll.map(
       lambda image: sampleGridCore(image, grid, scale)
    ).flatten()
  
    samples = samples.filter(ee.Filter.neq(control, None))
  
    if doexport:
        t = ee.batch.Export.table.toDrive(samples, fname, '')
        t.start()

        return {'samples':samples, 'task':t}

    else:
        return samples
