from gee_tools.seasonal import getS1Plus, getS2Plus
import ee


ee.Initialize()
# TODO: look into how to make this a unittest proper. Does CI work with EE?

# ROI around Nairobi, Kenya.
roi = ee.Geometry.Rectangle([
    36.7697414344982, -1.3481849306508715,
    37.095211405201326, -1.1278233090907424
])


s1seascoll = getS1Plus(roi, 2018, True, True)

s1seasimg = ee.Image(s1seascoll.first())

print "Test image for S1 has the following bands: "
print s1seasimg.bandNames().getInfo()
print '\n'


s2seascoll = getS2Plus(roi, 2018, True)

s2seasimg = ee.Image(s2seascoll.first())

print "Test image for S2 has the following bands: "
print s2seasimg.bandNames().getInfo()
print '\n'