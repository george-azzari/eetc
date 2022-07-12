"""
Author: George Azzari (gazzari@stanford.edu)
Center on Food Security and the Environment
Department of Earth System Science
Stanford University
"""

# pylint: disable=g-bad-import-order
# Using lowercase function naming to match the JavaScript names.
# pylint: disable=g-bad-name

import ee
import numpy as np
import gee_tools.datasources.optical_datasources as optix

from __future__ import print_function


def addTimeUnit(image, refdate):

    date = image.date()
  
    # Letting user pick the reference date
    dyear = date.difference(refdate, 'year')
    t = image.select(0).multiply(0).add(dyear).select([0],['t']).float()
    
    imageplus = image.addBands(t)

    return imageplus.set('timeunit', dyear)


def addConstant(image):

    return image.addBands(image.select(0).multiply(0).add(1).select([0], ['constant']))


def appendBand(current, previous):

    # Append it to the result (Note: only return current item on first element/iteration)
    accum = ee.Algorithms.If(ee.Algorithms.IsEqual(previous, None), current, current.addBands(ee.Image(previous)))
    # Return the accumulation
    return accum


def getNthHarmonic(n, omega, timeband):
    
    timerad = timeband.multiply(ee.Number(n)).multiply(ee.Number(2*np.pi)).multiply(ee.Number(omega))
    
    hcos = timerad.cos().rename(ee.String('cos').cat(ee.String(ee.Number(n).toInt())))
    hsin = timerad.sin().rename(ee.String('sin').cat(ee.String(ee.Number(n).toInt())))

    return ee.List([hcos, hsin])


def addHarmonics(image, timeband, omega, nharmonics):

    timeradians = ee.List.sequence(1, nharmonics, 1).map(lambda n: getNthHarmonic(n, omega, timeband)).flatten()

    # Convert list into a collection and smash into an image.
    timeradsimg = ee.Image(ee.ImageCollection.fromImages(timeradians).iterate(appendBand))
  
    constant = image.select(timeband).divide(image.select(timeband)).rename('constant')

    # Updating the input image w/ new harmonic bands
    image = image.addBands(timeradsimg).addBands(constant)
    image = image.set('independents', timeradsimg.bandNames().add(timeband).add('constant'))
  
    return image


def getHarmonicCollection(collection, omega, nharmonics, timeband, refdate):
    """
    Attach harmonic independents to an input image collection
    """

    if nharmonics is None:
        nharmonics = 2
    
    # If there is no time band in the collection generate the standard one
    if timeband is None:
        # Without reference date, use Jan-1 of the collection's first year.
        if refdate is None:
            d = ee.Image(collection.first()).date()
            refdate = ee.Date.fromYMD(d.get('year'), 1, 1)

        timeband = 't'
        collection = collection.map(lambda img: addTimeUnit(img, refdate)) 

    # Add harmonic terms as new image bands to collection.
    harmonic_coll = collection.map(lambda img: addHarmonics(img, timeband, omega, nharmonics))

    return  harmonic_coll


# def getDummyTime(startdate, enddate, ndays, dummyimg, addharmonics, omega, nharmonics):

#     diff = enddate.difference(startdate, 'day')
#     increments = ee.List.sequence(1, diff, null, ndays)
  
#     # Had to add the "round" part to avoid sub-daily increments when ndays is not 365.
#     dates = increments.map(lambda i: startdate.advance(ee.Number(i).round(), 'day'))


def arrayimgHarmonicRegr(harmonicoll, dependent, independents):

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


def imageHarmonicRegr(harmonicoll, dependent, independents):

    independents = ee.List(independents)
    dependent = ee.String(dependent)

    hregr = arrayimgHarmonicRegr(harmonicoll, dependent, independents)
    
    # New names for coefficients
    newnames = independents.map(lambda b: dependent.cat(ee.String('_')).cat(ee.String(b)))
    # Turn the array image into a multi-band image of coefficients.
    imgcoeffs = hregr.select('coefficients').arrayProject([0]).arrayFlatten([independents])
    imgcoeffs = imgcoeffs.select(independents, newnames)

    # The band 'residuals' the *root mean square* of the residuals (RMSE)
    rmse = hregr.select('residuals').arrayProject([0]).arrayFlatten([[dependent.cat(ee.String('_rmse'))]])

    # Computing variance and mean across collection (needed for R2)
    totreducer = ee.Reducer.sampleVariance().combine(ee.Reducer.count(), None, True)
    totreducer = totreducer.combine(ee.Reducer.mean(), None, True)
    stats = harmonicoll.select(dependent).reduce(totreducer)

    varbname = dependent.cat(ee.String('_variance'))
    variance = stats.select(varbname)
    
    # UPDATE: computing R2 from RMSE and variance can be done on the fly.
    #         Suppressing it here can save storage if asset is exported.
    # # Computing R2
    # r2bandn = dependent.cat(ee.String('_r2'))
    # r2 = ee.Image(1).updateMask(variance).subtract(rmse.pow(2).divide(variance)).select([0], [r2bandn])

    imgcoeffs = imgcoeffs.addBands(rmse).addBands(variance)#.addBands(r2)

    return imgcoeffs


def getHarmonicCoeffs(harmonicoll, bands, independents, ascoll):

    if independents is None:
        independents = ee.List(ee.Image(harmonicoll.first()).get('independents'))

    coefficients = ee.List(bands).map(lambda band: imageHarmonicRegr(harmonicoll, band, independents))
    coeffcoll = ee.ImageCollection.fromImages(coefficients)

    if ascoll:
        return coeffcoll
    else:
        return ee.Image(coeffcoll.iterate(appendBand))



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
        return ee.Image(coeffcoll.iterate(appendBand))


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
    lx = optix.LandsatSR(region, start_date, end_date).mergedqam
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
