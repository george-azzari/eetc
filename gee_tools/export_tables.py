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

    return ee.Feature(ee.Feature(None, s).copyProperties(image, None, image.bandNames()).copyProperties(point, None, image.bandNames())).set({

        #'PTLON': ee.Number(ee.List(point.geometry().coordinates()).get(0)),
        #'PTLAT': ee.Number(ee.List(point.geometry().coordinates()).get(1)),
        #'OBSDATE': ee.Date(image.get('system:time_start')).format(),
        #'MSTIME': image.get('system:time_start')

    })


def reducegrid_core(image, grid, scale):
    depth = image.bandNames().length()
    samples = grid.map(lambda point: get_point(image, point, scale, depth))
    return samples


def reducegrid_image(image, grid, scale, control, doexport, fname):

    samples = reducegrid_core(image, grid, scale)

    if doexport:
        t = export_features(samples, fname, export_to='drive')
        return {'samples':samples, 'task':t}

    else:
        return samples


def reducegrid_imgcoll(imagecoll, grid, scale, control, doexport, fname):

    samples = imagecoll.map(
       lambda image: reducegrid_core(image, grid, scale)
    ).flatten()
  
    samples = samples.filter(ee.Filter.neq(control, None))
  
    if doexport:
        t = export_features(samples, fname, export_to='drive')
        return {'samples':samples, 'task':t}

    else:
        return samples


def sampleregions_auto_image(image, regions, scale, controlvar, doexport, fname):

    samples = image.sampleRegions(
        collection=regions,
        properties=None,
        scale=scale,
        projection=None,
        tileScale=16)

    samples = samples.filter(ee.Filter.neq(controlvar, None))

    if doexport:
        task = export_features(samples, fname, export_to='drive')
        return {'samples': samples, 'task': task}

    else:
        return samples


def sampleregion_image(image, region, scale, npx):

    samples = image.sample(
        region=region.geometry(),
        scale=scale,
        projection=None,
        factor=None,
        numPixels=npx,
        seed=12345,
        dropNulls=True,
        tileScale=16)

    samples = samples.map(lambda p: ee.Feature(p.copyProperties(region).copyProperties(image)))

    return samples


def sampleregions_image(image, regions, scale, npx, doexport, fname):

    samples = regions.map(lambda region: sampleregion_image(image, region, scale, npx)).flatten()

    if doexport:
        task = export_features(samples, fname, export_to='drive')
        return {'samples': samples, 'task': task}

    else:
        return samples