from gee_tools.seasonal import getS1Plus
import ee


ee.Initialize()

# ROI around Nairobi, Kenya.
roi = ee.Geometry.Rectangle([
    36.7697414344982, -1.3481849306508715,
    37.095211405201326, -1.1278233090907424
])

# Generate seasonal collection for year 2018.
s1coll = getS1Plus(roi, 2018, correctlia=False, addspeckle=True)

testimg = ee.Image(s1coll.first())

print "Test image has the following bands: "
print testimg.bandNames().getInfo()