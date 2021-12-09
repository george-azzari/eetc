"""
Author: Shruti Jain

Obtains correction coefficients for converting Sentinel 2 TOA to SR
"""

import gee_tools.datasources.sentinel2 as s2_1c
import gee_tools.datasources.sentinel2_2a as s2_2a
from gee_tools.export_tables import export_features
import ee

ee.Initialize()

def _mergejoin(joinedelement):
    #  Inner join returns a FeatureCollection with a primary and secondary set of
    #  properties. Properties are collapsed into different bands of an image.
    return ee.Image.cat(joinedelement.get('primary'), joinedelement.get('secondary'))


#  Convenience function for joining two collections based on system:time_start
def joincoll(coll1, coll2):
    eqfilter = ee.Filter.equals(rightField='system:index', leftField='system:index')
    join = ee.Join.inner()
    joined = ee.ImageCollection(join.apply(coll1, coll2, eqfilter))
    return joined.map(_mergejoin).sort('system:index')


def get_samples(img, num_pixels, geometry, mask):
    
    img = img.clip(geometry)
    img = img.updateMask(mask)
    points = img.sample(
        scale=10, 
        numPixels=num_pixels,
        seed=23,
        geometries=True
        )
    
    return points

def get_coeffs(sample_points, bandname, country):
    coeffs = sample_points.reduceColumns(
        reducer=ee.Reducer.linearFit(), 
        selectors=[bandname+'1', bandname+'2']
        )
    n = sample_points.size()
    coeffs = ee.Feature(None, coeffs).set('n', n, 'band', bandname,'country', country)
    return coeffs

# For s2 toa
def maskClouds_toa(img, bandnames):
    qa60 = img.select('QA60_DECODED')
    clear = qa60.updateMask(qa60.eq(1))
    img = img.updateMask(clear)
    bandnames_new = [b+'1' for b in bandnames]
    return img.select(bandnames, bandnames_new)

# For s2 sr
def maskClouds_sr(img, bandnames):
    scl = img.select(['SCL'])
    clear = scl.updateMask(scl.eq(4).Or(scl.eq(5)))
    img = img.updateMask(clear)
    bandnames_new = [b+'2' for b in bandnames]
    return img.select(bandnames, bandnames_new)

def get_corr_coeffs(geometry, country, max_slope=10, start_date=ee.Date('2019-01-01'), end_date=ee.Date('2019-12-31'), num_pixels=2):
    srtm = ee.Image("USGS/SRTMGL1_003")
    topog = ee.Algorithms.Terrain(srtm).select(['elevation', 'slope', 'aspect'],['ELEV', 'SLO', 'ASP'])
    slo = topog.select('SLO')
    mask = slo.lt(max_slope)
    mask = mask.clip(geometry)

    bandnames = ['AEROS','BLUE','GREEN','RED','RDED1','RDED2','RDED3','NIR','RDED4','VAPOR','SWIR1','SWIR2']
    s2coll_1c = s2_1c.Sentinel2TOA(geometry, start_date, end_date, addVIs=False, addCloudMasks=True).get_img_coll()
    s2coll_1c = s2coll_1c.map(lambda img: maskClouds_toa(img, bandnames))
    s2coll_2a = s2_2a.Sentinel2SR(geometry, start_date, end_date, addVIs=False, addCloudMasks=True).get_img_coll()
    s2coll_2a = s2coll_2a.map(lambda img: maskClouds_sr(img, bandnames))
    imgcoll = joincoll(s2coll_1c, s2coll_2a)
    imgcoll = ee.ImageCollection(imgcoll.toList(imgcoll.size()))

    sample_points = ee.FeatureCollection(imgcoll.map(lambda img: get_samples(img, num_pixels, geometry, mask))).flatten()
    coefficients = ee.FeatureCollection(ee.List([
        get_coeffs(sample_points, bandname, country) for bandname in bandnames
    ]))

    task = export_features(coefficients, f"s2_correction_coeffs/s2_corr_coeffs_{country}.csv")
    
    return task

    