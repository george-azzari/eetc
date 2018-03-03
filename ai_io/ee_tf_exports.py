import ee


def _sample_patch(point, patchesarray, scale):
    arrays_samples = patchesarray.sample(
        region=point.geometry(),
        scale=scale,
        #           projection='EPSG:32610',
        projection='EPSG:3857',
        factor=None,
        numPixels=None,
        dropNulls=False,
        tileScale=12

    )

    #     return ee.Feature(arrays_samples.copyProperties(point))
    return arrays_samples


def export_patches(img, scale, ksize, points, doexport, tocloud, selectors, dropselectors, mybucket, prefix, fname):

    kern = ee.Kernel.square(ksize, 'pixels')
    patches_array = img.neighborhoodToArray(kern)

    # sampleRegions does not cut it for larger collections; using mapped sample instead.
    patches_samps = points.map(lambda pt: _sample_patch(pt, patches_array, scale)).flatten()

    if doexport:
        # Export to a TFRecord file in Cloud Storage, creating a file
        # at gs://mybucket/prefix/fname.tfrecord
        # which you can load directly in TensorFlow.

        if selectors is None:
            selectors = ee.Feature(patches_samps.first()).propertyNames()

        if dropselectors is not None:
            selectors = selectors.removeAll(dropselectors)

        if tocloud:
            task = ee.batch.Export.table.toCloudStorage(

                collection=patches_samps,
                description=fname,
                bucket=mybucket,
                fileNamePrefix=prefix + fname,
                #                 fileFormat='CSV',
                fileFormat='TFRecord',
                selectors=selectors,
            )

        else:
            task = ee.batch.Export.table.toDrive(

                collection=patches_samps,
                description=fname,
                folder='',
                fileNamePrefix=None,
                # fileFormat= 'CSV',
                fileFormat='TFRecord',
                selectors=selectors

            )

        task.start()

    return patches_samps