"""
Author: George Azzari (gazzari@stanford.edu)
Center on Food Security and the Environment
Department of Earth System Science
Stanford University
"""

import ee
import numpy as np


def append_band(current, previous):

    # Append it to the result (Note: only return current item on first element/iteration)
    accum = ee.Algorithms.If(ee.Algorithms.IsEqual(previous, None), current, current.addBands(ee.Image(previous)))
    # Return the accumulation
    return accum


def add_timevars(image, tunit):

    date = image.date()
    doy = date.getRelative('day', 'year').add(1)
    doyimage = image.select(0).multiply(0).add(doy).select([0], ['DOY']).toInt16()

    month = date.getRelative('month', 'year')
    monthimage = image.select(0).multiply(0).add(month).select([0], ['MONTH']).toInt16()

    dyear = date.difference(ee.Date('1970-01-01'), 'year')
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

    timeRadians = image.select('t').multiply(2 * np.pi * omega)

    cost = timeRadians.cos().rename(['cos'])
    sint = timeRadians.sin().rename(['sin'])
    sincost = timeRadians.sin().multiply(timeRadians.cos()).rename(['sincos'])

    return image.addBands(cost).addBands(sint).addBands(sincost)


def get_harmonic_coll(collection, omega):

    omega = 1.6

    f = lambda img: add_harmonics(img, omega)

    # Add harmonic terms as new image bands.
    harmonic_coll = collection.map(f)

    return harmonic_coll


def arrayimg_harmon_regr(harmonicoll, dependent, independents):


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

    regressors= harmonicoll.select(independents.add(dependent))
    regression = regressors.reduce(ee.Reducer.linearRegression(independents.length(), 1))

    return  regression


def get_prediction(harmonicimg, regrcoeffimg, dependent, independents):

    harmonicimg = ee.Image(harmonicimg).select(independents)

    predicted = harmonicimg.multiply(regrcoeffimg).reduce(ee.Reducer.sum())
    predicted =predicted.select([0], [ee.String('pred_').cat(dependent)])

    return predicted


def get_residual(harmonicimg, regrcoeffimg, dependent, independents):

    predicted = get_prediction(harmonicimg, regrcoeffimg, dependent, independents)
    residual = harmonicimg.select(dependent).subtract(predicted).pow(2)

    return residual


def image_harmon_regr(harmonicoll, dependent, independents, myrmse=False):

    hregr = arrayimg_harmon_regr(harmonicoll, dependent, independents)

    independents = ee.List(independents)
    dependent = ee.String(dependent)

    totreducer = ee.Reducer.sampleVariance()
    totreducer = totreducer.combine(ee.Reducer.count(), None, True)
    totreducer = totreducer.combine(ee.Reducer.mean(), None, True)

    stats = harmonicoll.select(dependent).reduce(totreducer)

    # New names for coefficients
    newnames = independents.map(lambda b: dependent.cat(ee.String('_')).cat(ee.String(b)))

    # Turn the array image into a multi-band image of coefficients.
    imgcoeffs = hregr.select('coefficients').arrayProject([0]).arrayFlatten([independents])
    imgcoeffs = imgcoeffs.select(independents, newnames)

    if myrmse:

        rss = harmonicoll.map(lambda himg: get_residual(himg, imgcoeffs, dependent, independents)).sum()
        rss = rss.select([0],[dependent.cat(ee.String('_rss'))])

        count = stats.select(dependent.cat(ee.String('_count')))
        rmse = rss.divide(count).sqrt().select([0], [dependent.cat(ee.String('_rmse2'))])

    else:
        # The band 'residuals' the *root mean square* of the residuals (RMSE)
        rmse = hregr.select('residuals').arrayProject([0]).arrayFlatten([[dependent.cat(ee.String('_rmse'))]])

    variance = stats.select(dependent.cat(ee.String('_variance')))

    r2bandn =dependent.cat(ee.String('_r2'))
    r2 = ee.Image(1).updateMask(variance).subtract(rmse.pow(2).divide(variance)).select([0], [r2bandn])

    return imgcoeffs.addBands(stats).addBands(rmse).addBands(r2)


def allbands_harmon_regr(collection, bands, independents, ascoll):

    # Add harmonic terms as new image bands.
    hcoll = get_harmonic_coll(collection, 1.6)

    coefficients = ee.List(bands).map(lambda dependent: image_harmon_regr(hcoll, dependent, independents))
    coeffcoll = ee.ImageCollection.fromImages(coefficients)

    if ascoll:
        return coeffcoll

    else:
        return ee.Image(coeffcoll.iterate(append_band))



