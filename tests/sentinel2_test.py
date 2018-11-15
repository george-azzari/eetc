import unittest

import ee

from gee_tools.datasources import sentinel2


ee.Initialize()
# TODO: look into how to make this a unittest proper. Does CI work with EE?

# ROI around Nairobi, Kenya.
roi = ee.Geometry.Rectangle([
    36.7697414344982, -1.3481849306508715,
    37.095211405201326, -1.1278233090907424
])


# Define reference dates.
sdate = ee.Date('2017-11-1')
edate = ee.Date('2018-11-1')


"""
Test base collection (raw S1 assets).
"""
# Generate base collection.
s2base = sentinel2.Sentinel2TOA(roi, sdate, edate)
s2base.build_img_coll()

testimg_base = ee.Image(s2base.coll.first())

print "Test image base has the following bands: "
print testimg_base.bandNames().getInfo()
print '\n'


"""
Test VIs.
"""
# Get collection with extra bands.
s2extra = s2base.get_img_coll(addVIs=True, addCloudMasks=False)

testimg_s2extra = ee.Image(s2extra.first())

print "Test image with VIs has the following bands: "
print testimg_s2extra.bandNames().getInfo()
print '\n'


"""
Test Cloud masks.
"""
# Get collection with extra bands.
s2qa = s2base.get_img_coll(addVIs=True, addCloudMasks=True)

testimg_s2qa = ee.Image(s2qa.first())

print "Test image with QA  has the following bands: "
print testimg_s2qa.bandNames().getInfo()
print '\n'
