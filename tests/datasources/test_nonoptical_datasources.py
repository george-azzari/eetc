"""
"""
import unittest

import ee

from gee_tools.datasources import nonoptical_datasources as nonop
from tests.test_utils import compare_bands


class ModisDailyLstLongRunningAverageTestCase(unittest.TestCase):

    def setUp(self):
        ee.Initialize()

        # ROI around Nairobi, Kenya.
        self.roi = ee.Geometry.Rectangle([
            36.7697414344982, -1.3481849306508715,
            37.095211405201326, -1.1278233090907424
        ])

        # Define reference dates.
        self.sdate = ee.Date('2017-11-1')
        self.edate = ee.Date('2018-11-1')

    def test_basic(self):
        season_to_int_months = {'WINTER': [12, 1, 2]}
        datasource = nonop.ModisDailyLstLongRunningAverage(season_to_int_months=season_to_int_months)
        img = ee.Image(datasource.get_img_coll().first())
        band_names = sorted(img.bandNames().getInfo())
        expected_band_names = nonop.ModisDailyLstLongRunningAverage.band_names(season_to_int_months=season_to_int_months)
        expected_band_names = sorted(expected_band_names)
        self.assertListEqual(band_names, expected_band_names)


if __name__ == '__main__':
    unittest.main()
