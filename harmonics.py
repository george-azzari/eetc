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
    dyimage = image.select(0).multiply(0).add(dyear).select([0] ,['DYEAR']).float()

    ms = ee.Number(image.get('system:time_start'))
    msimage = image.select(0).multiply(0).add(ms).select([0] ,['MSTIME'])

    imageplus = ee.Image.cat([image, doyimage, monthimage, msimage, dyimage])
    imageplus = imageplus.set({'DOY': doy, 'MONTH': month, 'YEAR': date.get('year'), 'DYEAR':dyear})

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


def get_coeffs(band, harmonicLandsat, harmonicIndependents):

    # Name of the dependent variable.
    dependent = ee.String(band)
    # The output of the regression reduction is a 4x1 array image.
    regression_coeff = harmonicLandsat.select(harmonicIndependents.add(dependent)).reduce(
        ee.Reducer.linearRegression(harmonicIndependents.length(), 1))

    # Turn the array image into a multi-band image of coefficients.
    regression_imgcoeff = regression_coeff.select('coefficients').arrayProject([0]).arrayFlatten(
        [harmonicIndependents]).select(harmonicIndependents,
                                       harmonicIndependents.map(
                                           lambda b: dependent.cat(ee.String('_')).cat(ee.String(b))))

    return regression_imgcoeff


def get_harmonic_coeffs(collection, bands, ascoll):

    # Use these independent variables in the harmonic regression.
    harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin', 'sincos'])

    # Add harmonic terms as new image bands.
    harmonicLandsat = get_harmonic_coll(collection, omega=1.6)

    f = lambda band: get_coeffs(band, harmonicLandsat, harmonicIndependents)
    coefficients = ee.List(bands).map(f)

    coeffcoll = ee.ImageCollection.fromImages(coefficients)

    if ascoll:
        return coeffcoll
    else:
        return ee.Image(coeffcoll.iterate(append_band))