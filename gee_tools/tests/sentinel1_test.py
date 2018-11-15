from gee_tools.datasources import sentinel1
import ee


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
s1base = sentinel1.Sentinel1(roi, sdate, edate)
s1base.build_img_coll()

testimg_base = ee.Image(s1base.coll.first())

print "Test image base has the following bands: "
print testimg_base.bandNames().getInfo()
print '\n'


"""
Test extra bands (difference and polarization).
"""
# Get collection with extra bands.
s1extra = s1base.get_img_coll(correctlia=False, addbands=True, addspeckle=False, addtexture=False, orbit='ascending')

testimg_s1extra = ee.Image(s1extra.first())

print "Test image with extras has the following bands: "
print testimg_s1extra.bandNames().getInfo()
print '\n'


"""
Test speckle-correction.
"""
# Get collection with extra bands.
s1speckle = s1base.get_img_coll(correctlia=False, addbands=True, addspeckle=True, addtexture=False, orbit='ascending')

testimg_s1speckle = ee.Image(s1speckle.first())

print "Test image with speckle-correction has the following bands: "
print testimg_s1speckle.bandNames().getInfo()
print '\n'


"""
Test GLCM texture.
"""
# Get collection with extra bands.
s1glcm = s1base.get_img_coll(correctlia=False, addbands=True, addspeckle=True, addtexture=True, orbit='ascending')

testimg_s1glcm = ee.Image(s1glcm.first())

print "Test image with speckle-correction and GLCM texture has the following bands: "
print testimg_s1glcm.bandNames().getInfo()
print '\n'



