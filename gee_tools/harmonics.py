"""
Author: George Azzari (gazzari@stanford.edu)
Center on Food Security and the Environment
Department of Earth System Science
Stanford University
"""

import ee
import numpy as np
import gee_tools.datasources.optical_datasources as optix

def add_timeunit(image, refdate):
    
    date = image.date()
    
    #Letting user pick the reference date
    dyear = date.difference(refdate, 'year')
    t = image.select(0).multiply(0).add(dyear).select([0],['t']).float()
    
    imageplus = image.addBands(t)
    return imageplus.set('timeunit', dyear)

def append_band(current, previous):

    # Append it to the result (Note: only return current item on first element/iteration)
    accum = ee.Algorithms.If(ee.Algorithms.IsEqual(previous, None), current, current.addBands(ee.Image(previous)))
    # Return the accumulation
    return accum

def add_constant(image):

    return image.addBands(image.select(0).multiply(0).add(1).select([0], ['constant']))

def add_harmonics(image, timeband, omega, nharmonics=2):
    
    def _add_harmonic(n):
        timerad = image.select(timeband).multiply(n * 2 * np.pi * omega)
        cos = timerad.cos().rename([f"cos{n}"])
        sin = timerad.sin().rename([f"sin{n}"])
        return ee.List([cos, sin])
    
    timeradians = ee.List([_add_harmonic(n) for n in range(1, nharmonics+1)])
    timeradians = timeradians.flatten()
    
    # Convert list into a collection and smash into an image.
    timeradsimg = ee.Image(ee.ImageCollection.fromImages(timeradians).iterate(append_band))
    constant = image.addBands(image.select(0).multiply(0).add(1).select([0], ['constant']))
    
    return image.addBands(timeradsimg).addBands(constant).set('independents', timeradsimg.bandNames().add(timeband).add('constant'))

def get_harmonic_coll(collection, omega=1, nharmonics=2, timeband=None, refdate=None):
    
    if timeband is None:
        if refdate is None:
            d = ee.Image(collection.first()).date()
            refdate = ee.Date.fromYMD(d.get('year'), 1, 1)
        timeband = 't'
        collection = collection.map(lambda img: add_timeunit(img, refdate))
    
    # Add harmonic terms as new image bands.
    f = lambda img: add_harmonics(img, timeband, omega, nharmonics)
    harmonic_coll = collection.map(f)

    return harmonic_coll


""" Generate a dummy collection of harmonis indepents with a cadence of ndays.
    Useful to generate "predicted" values from pre-fitted coefficients """
def get_dummy_time(startdate, enddate, ndays, dummyimg, addharmonics, omega=1, nharmonics=2):
    
    def _add_dummy(d):
        d = ee.Date(d)
        dummy = dummyimg.add(d.millis()).select([0],['dummy']).set('date', d.format(), 'system:time_start', d.millis())
        dummy = add_timeunit(dummy, startdate)
        return dummy

    diff = enddate.difference(startdate, 'day')
    increments = ee.List.sequence(1, diff, None, ndays)
    
    #Had to add the "round" part to avoid sub-daily increments when ndays is not 365.
    dates = increments.map(lambda i: startdate.advance(ee.Number(i).round(), 'day'))

    dummys = dates.map(_add_dummy)
    
    imgcoll =  ee.ImageCollection.fromImages(dummys).set('startdate', startdate, 'enddate', enddate, 'ndays', ndays)
    
    if addharmonics:
        imgcoll = get_harmonic_coll(imgcoll, omega, nharmonics, 't')  
    
    return  imgcoll

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

    return regression.reproject(ee.Image(harmonicoll.first()).projection())


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
    totreducer = ee.Reducer.sampleVariance().combine(ee.Reducer.count(), None, True).combine(ee.Reducer.mean(), None, True)
    stats = harmonicoll.select(dependent).reduce(totreducer)
    stats = stats.reproject(ee.Image(harmonicoll.first()).projection())
    variance = stats.select(dependent.cat(ee.String('_variance')))

    # Computing R2
    r2bandn = dependent.cat(ee.String('_r2'))
    r2 = ee.Image(1).updateMask(variance).subtract(rmse.pow(2).divide(variance)).select([0], [r2bandn])
    
    return imgcoeffs.addBands(stats).addBands(rmse).addBands(r2)


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

    if not independents:
        independents = ee.List(ee.Image(harmonicoll.first()).get('independents'))

    f= lambda dependent: hrmregr_single(harmonicoll, dependent, independents)
    coeffcoll = ee.ImageCollection.fromImages(ee.List(dependents).map(f))

    if ascoll:
        return coeffcoll

    else:
        return ee.Image(coeffcoll.iterate(append_band))

def run_std_regressions(imagecoll, bands, refdate, independents=None, omega=1, nharmonics=2, ascoll=False):
    
    # Init a new time variable with default reference date 
    harmcoll = get_harmonic_coll(imagecoll, omega, nharmonics, None, refdate)
    
    # Compute Regression Coefficients
    harmcoeffs = hrmregr_multi(harmcoll, bands, independents, ascoll)
    
    return  harmcoeffs

def fit_harmonics(harmcoeffs, imgcoll, omega, nharmonics, bands, refdate):

    def _fit_harmonic(image):
        indep = ee.List(ee.Image(image).get('independents'))
        def _fit_band(band):
            indepcoeffs = indep.map(lambda s: ee.String(band).cat(ee.String("_")).cat(s))
            return image.select(indep).multiply(harmcoeffs.select(indepcoeffs)).reduce('sum').rename(ee.String(band).cat(ee.String("_HARMFIT")))
        fittedimg = ee.List(bands).map(_fit_band)
        fittedimg = ee.ImageCollection.fromImages(fittedimg.add(image)).iterate(append_band)
        return fittedimg

    harmcoll = get_harmonic_coll(imgcoll, omega, nharmonics, None, refdate)

    # Compute fitted values.
    fittedcoll = harmcoll.map(_fit_harmonic)
    return fittedcoll

# NOTE: Out of date
# def lx_hregr(region, start_date, end_date, omega=1.5, imgmask=None, bands=None, rmbands=None,
#              independents=None, addcount=True):
#     """
#     Generate harmonics composite for a merged Landsat SR collection.
#     :param region: ee.Feature
#     :param start_date: ee.Date
#     :param end_date: ee.Date
#     :param omega: the omega factor for the Fourier series
#     :param imgmask: mask to use on individual images prior generating the composite (optional)
#     :param bands: the bands to use as dependent variables (optional, defaults to all optical bands and indexes)
#     :param rmbands: bands to remove from default (optional)
#     :param addcount: whether a count band should be added
#     :return: ee.Image (composite)
#     """
#     if independents is None:
#         # NOTE: removed 't' (linear term)
#         independents = ee.List(['constant', 'cos', 'sin', 'cos2', 'sin2'])

#     # TODO: update to new Collection1 collections
#     lx = optix.LandsatSR(region, start_date, end_date).mergedqam
#     lx = lx.select(['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']).map(optix.addVIs)

#     if imgmask is not None:
#         lx = lx.map(lambda img: img.updateMask(imgmask))

#     hlx = get_harmonic_coll(lx, omega)

#     if bands is None:
#         nonoptical = ee.List(['t', 'DOY', 'MONTH', 'MSTIME', 'DYEAR', 'constant'])
#         bands = ee.Image(lx.first()).bandNames().removeAll(nonoptical)

#         if rmbands is not None:
#             bands = bands.removeAll(ee.List(rmbands))

#     else:
#         bands = ee.List(bands)

#     allcoeffs = hrmregr_multi(hlx, bands, independents, False)
#     allcoeffs = allcoeffs.set('omega', omega,
#                               'n', 2,
#                               'start_date', start_date.millis(),
#                               'end_date', end_date.millis(),
#                               'formula', 'A + Bt + Ccos(2piwt) + Dsin(2piwt) + Ecos(4piwt) + Fsin(4piwt)')

#     if addcount:
#         count = lx.select(['NIR']).count().select([0], ['count'])
#         allcoeffs = allcoeffs.addBands(count)

#     return allcoeffs
