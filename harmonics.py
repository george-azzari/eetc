"""
Author: George Azzari (gazzari@stanford.edu)
Center on Food Security and the Environment
Department of Earth System Science
Stanford University
"""

import ee
import numpy as np
import optical_datasources as optix


def append_band(current, previous):

    # Append it to the result (Note: only return current item on first element/iteration)
    accum = ee.Algorithms.If(ee.Algorithms.IsEqual(previous, None), current, current.addBands(ee.Image(previous)))
    # Return the accumulation
    return accum


def add_timevars(image, tunit):

    # TODO: update this to be a simpler add_timeunit. Other info are irrelevant to the regression.

    date = image.date()
    doy = date.getRelative('day', 'year').add(1)
    doyimage = image.select(0).multiply(0).add(doy).select([0], ['DOY']).toInt16()

    month = date.getRelative('month', 'year')
    monthimage = image.select(0).multiply(0).add(month).select([0], ['MONTH']).toInt16()

    jan1 = ee.Date.fromYMD(date.get('year'), 1, 1)
    dyear = date.difference(jan1, 'year')
    dyimage = image.select(0).multiply(0).add(dyear).select([0], ['DYEAR']).float()

    ms = ee.Number(image.get('system:time_start'))
    msimage = image.select(0).multiply(0).add(ms).select([0], ['MSTIME'])

    imageplus = ee.Image.cat([image, doyimage, monthimage, msimage, dyimage])
    imageplus = imageplus.set({'DOY': doy, 'MONTH': month, 'YEAR': date.get('year'), 'DYEAR': dyear})

    bnames = ee.List(imageplus.bandNames()).replace(tunit, 't')

    return imageplus.rename(bnames)


def add_constant(image):

    return image.addBands(image.select(0).multiply(0).add(1).select([0], ['constant']))


def add_harmonics(image, omega):

    image = add_timevars(image, 'DYEAR')
    image = add_constant(image)

    timeRadians = image.select('t').multiply(2 * np.pi * omega)
    timeRadians2 = image.select('t').multiply(4 * np.pi * omega)

    cost = timeRadians.cos().rename(['cos'])
    sint = timeRadians.sin().rename(['sin'])
    cost2 = timeRadians2.cos().rename(['cos2'])
    sint2 = timeRadians2.sin().rename(['sin2'])

    return image.addBands(cost).addBands(sint).addBands(cost2).addBands(sint2)


def get_harmonic_coll(collection, omega=1.5):

    f = lambda img: add_harmonics(img, omega)

    # Add harmonic terms as new image bands.
    harmonic_coll = collection.map(f)

    return harmonic_coll


def _arrayimg_hrmregr_single(harmonicoll, dependent, independents):

    """
     The first output is a coefficients array with dimensions (numX, numY)
         each column contains the coefficients for the corresponding dependent variable.
         The second output is a vector of the *root mean square* of the residuals of each dependent variable.
         Both outputs are null if the system is underdetermined, e.g. the number of inputs is less than or equal to numX.

     :param harmonicoll:
     :param dependent:
     :param independents:
     :return:
     """
    independents = ee.List(independents)
    dependent = ee.String(dependent)

    regressors = harmonicoll.select(independents.add(dependent))
    regression = regressors.reduce(ee.Reducer.linearRegression(independents.length(), 1))

    return regression


def hrmregr_single(harmonicoll, dependent, independents):

    independents = ee.List(independents)
    dependent = ee.String(dependent)

    hregr = _arrayimg_hrmregr_single(harmonicoll, dependent, independents)

    # New names for coefficients
    newnames = independents.map(lambda b: dependent.cat(ee.String('_')).cat(ee.String(b)))
    # Turn the array image into a multi-band image of coefficients.
    imgcoeffs = hregr.select('coefficients').arrayProject([0]).arrayFlatten([independents])
    imgcoeffs = imgcoeffs.select(independents, newnames)

    # The band 'residuals' the *root mean square* of the residuals (RMSE)
    rmse = hregr.select('residuals').arrayProject([0]).arrayFlatten([[dependent.cat(ee.String('_rmse'))]])

    # Computing variance for R2
    totreducer = ee.Reducer.sampleVariance()
    variance = harmonicoll.select(dependent).reduce(totreducer)

    # UPDATE: computing R2 from RMSE and variance can be done on the fly.
    #         Suppressing it here can save storage if asset is exported.
    # # Computing R2
    # r2bandn = dependent.cat(ee.String('_r2'))
    # r2 = ee.Image(1).updateMask(variance).subtract(rmse.pow(2).divide(variance)).select([0], [r2bandn])

    imgcoeffs = imgcoeffs.addBands(rmse).addBands(variance)#.addBands(r2)

    return imgcoeffs


def _get_prediction(harmonicimg, regrcoeffimg, dependent, independents):

    harmonicimg = ee.Image(harmonicimg).select(independents)

    predicted = harmonicimg.multiply(regrcoeffimg).reduce(ee.Reducer.sum())
    predicted = predicted.select([0], [ee.String('pred_').cat(dependent)])

    return predicted


def _get_residual(harmonicimg, regrcoeffimg, dependent, independents):

    predicted = _get_prediction(harmonicimg, regrcoeffimg, dependent, independents)
    residual = harmonicimg.select(dependent).subtract(predicted).pow(2)

    return residual


def hrmregr_multi(harmonicoll, dependents, independents, ascoll):

    """
    :param harmonicoll: collection containing all bands of interest plus the independents.
    :param dependents: list of image bands on which regression is to be performed.
    :param independents: list of time bands (constant, time, sines, and cosines)
    :param ascoll: whether the result should be returned as a collection, rather than an image.
    :return: ee.Image containing coefficients for all dependents and independents.
    """

    f= lambda dependent: hrmregr_single(harmonicoll, dependent, independents)
    coeffcoll = ee.ImageCollection.fromImages(ee.List(dependents).map(f))

    if ascoll:
        return coeffcoll

    else:
        return ee.Image(coeffcoll.iterate(append_band))


def lx_hregr(region, start_date, end_date, omega=1.5, imgmask=None, bands=None, rmbands=None,
             independents=None, addcount=True):
    """
    Generate harmonics composite for a merged Landsat SR collection.
    :param region: ee.Feature
    :param start_date: ee.Date
    :param end_date: ee.Date
    :param omega: the omega factor for the Fourier series
    :param imgmask: mask to use on individual images prior generating the composite (optional)
    :param bands: the bands to use as dependent variables (optional, defaults to all optical bands and indexes)
    :param rmbands: bands to remove from default (optional)
    :param addcount: whether a count band should be added
    :return: ee.Image (composite)
    """
    if independents is None:
        # NOTE: removed 't' (linear term)
        independents = ee.List(['constant', 'cos', 'sin', 'cos2', 'sin2'])

    # TODO: update to new Collection1 collections
    lx = optix.LandsatSR(region, start_date, end_date).mergedcfm
    lx = lx.select(['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']).map(optix.addVIs)

    if imgmask is not None:
        lx = lx.map(lambda img: img.updateMask(imgmask))

    hlx = get_harmonic_coll(lx, omega)

    if bands is None:
        nonoptical = ee.List(['t', 'DOY', 'MONTH', 'MSTIME', 'DYEAR', 'constant'])
        bands = ee.Image(lx.first()).bandNames().removeAll(nonoptical)

        if rmbands is not None:
            bands = bands.removeAll(ee.List(rmbands))

    else:
        bands = ee.List(bands)

    allcoeffs = hrmregr_multi(hlx, bands, independents, False)
    allcoeffs = allcoeffs.set('omega', omega,
                              'n', 2,
                              'start_date', start_date.millis(),
                              'end_date', end_date.millis(),
                              'formula', 'A + Bt + Ccos(2piwt) + Dsin(2piwt) + Ecos(4piwt) + Fsin(4piwt)')

    if addcount:
        count = lx.select(['NIR']).count().select([0], ['count'])
        allcoeffs = allcoeffs.addBands(count)

    return allcoeffs
