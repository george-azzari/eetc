"""
# TODO: look into how to make this a unittest proper. Does CI work with EE?
"""
import unittest

import ee

from gee_tools.seasonal import getS1Plus, getS2Plus, sentinel_combined_medians


def compare_bands(testcase, img, expected_bands, assert_equal_kwargs=None):
    if assert_equal_kwargs is None:
        assert_equal_kwargs = {}

    img_bands = img.bandNames().getInfo()
    expected_bands = sorted(list(set(expected_bands)))
    img_bands = sorted(list(img_bands))

    testcase.assetEqual(img_bands, expected_bands, **assert_equal_kwargs)


class SeasonalTestCase(unittest.TestCase):

    def setUp(self):
        ee.Initialize()

        # ROI around Nairobi, Kenya.
        self.roi = ee.Geometry.Rectangle([
            36.7697414344982, -1.3481849306508715,
            37.095211405201326, -1.1278233090907424
        ])

    def test_seasonal(self):
        glcmvars = ['contrast', 'corr', 'var', 'savg', 'prom']
        seasmed = sentinel_combined_medians(self.roi, 2018, False, False, True, True, glcmvars)
        expected_bands = []  # TODO
        compare_bands(self, seasmed, expected_bands, {'msg': 'Test image for sentinel_combined_medians had the wrong bands'})

    def test_s1(self):
        s1seascoll = getS1Plus(self.roi, 2018, True, True)
        s1seasimg = ee.Image(s1seascoll.first())

        expected_bands = []  # TODO
        compare_bands(self, s1seasimg, expected_bands, {'msg': 'Test image for S1 had the wrong bands'})

    def test_s2(self):
        s2seascoll = getS2Plus(self.roi, 2018, True)
        s2seasimg = ee.Image(s2seascoll.first())

        expected_bands = []  # TODO
        compare_bands(self, s2seasimg, expected_bands, {'msg': 'Test image for S2 had the wrong bands'})


if __name__ == '__main__':
    unittest.main()
