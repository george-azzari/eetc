"""
# TODO: look into how to make this a unittest proper. Does CI work with EE?
"""
import unittest

import ee

from gee_tools.datasources import sentinel2
from tests.test_utils import compare_bands


class Sentinel2TestCase(unittest.TestCase):

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

    def test_base_collection(self):
        """
        Test base collection (raw S2 assets).
        """
        # Generate base collection.
        s2base = sentinel2.Sentinel2TOA(self.roi, self.sdate, self.edate)
        testimg_base = ee.Image(s2base.coll.first())  # TODO tests should be agnostic to the underlying implementation

        expected_bands = [
            u'AEROS', u'BLUE', u'GREEN', u'RED', u'RDED1', u'RDED2', u'RDED3',
            u'NIR', u'RDED4', u'VAPOR', u'CIRRU', u'SWIR1', u'SWIR2', u'QA10',
            u'QA20', u'QA60'
        ]
        compare_bands(self, testimg_base, expected_bands, {'msg': 'Sentinel 2 base image had the wrong bands'})

    def test_VIs(self):
        """
        Test VIs.
        """
        # Get collection with extra bands.
        s2base = sentinel2.Sentinel2TOA(self.roi, self.sdate, self.edate)
        s2extra = s2base.get_img_coll(addVIs=True, addCloudMasks=False)
        testimg_s2extra = ee.Image(s2extra.first())

        expected_bands = [
            u'AEROS', u'BLUE', u'GREEN', u'RED', u'RDED1', u'RDED2', u'RDED3', u'NIR', u'RDED4',
            u'VAPOR', u'CIRRU', u'SWIR1', u'SWIR2', u'QA10', u'QA20', u'QA60', u'NBR1', u'NBR2',
            u'STI', u'NDTI', u'CRC', u'REIP', u'GCVI', u'RDGCVI1', u'RDGCVI2', u'MTCI', u'MTCI2',
            u'WDRVI', u'GRWDRVI', u'RDWDRVI', u'RDNDVI1', u'RDNDVI2', u'NDVI'
        ]
        compare_bands(self, testimg_s2extra, expected_bands, {'msg': 'Sentinel 2 with VIs had the wrong bands'})

    def test_VIs_constructor(self):
        """
        Test VIs.
        """
        # Get collection with extra bands.
        s2base = sentinel2.Sentinel2TOA(
            self.roi, self.sdate, self.edate,
            addVIs=True, addCloudMasks=False,
        )
        s2extra = s2base.get_img_coll()
        testimg_s2extra = ee.Image(s2extra.first())

        expected_bands = [
            u'AEROS', u'BLUE', u'GREEN', u'RED', u'RDED1', u'RDED2', u'RDED3', u'NIR', u'RDED4',
            u'VAPOR', u'CIRRU', u'SWIR1', u'SWIR2', u'QA10', u'QA20', u'QA60', u'NBR1', u'NBR2',
            u'STI', u'NDTI', u'CRC', u'REIP', u'GCVI', u'RDGCVI1', u'RDGCVI2', u'MTCI', u'MTCI2',
            u'WDRVI', u'GRWDRVI', u'RDWDRVI', u'RDNDVI1', u'RDNDVI2', u'NDVI'
        ]
        compare_bands(self, testimg_s2extra, expected_bands, {'msg': 'Sentinel 2 with VIs had the wrong bands'})

    def test_cloud_mask(self):
        """
        Test Cloud masks.
        """
        # Get collection with extra bands.
        s2base = sentinel2.Sentinel2TOA(self.roi, self.sdate, self.edate)
        s2qa = s2base.get_img_coll(addVIs=True, addCloudMasks=True)
        testimg_s2qa = ee.Image(s2qa.first())

        expected_bands = [
            u'AEROS', u'BLUE', u'GREEN', u'RED', u'RDED1', u'RDED2', u'RDED3', u'NIR', u'RDED4', u'VAPOR',
            u'CIRRU', u'SWIR1', u'SWIR2', u'QA10', u'QA20', u'QA60', u'QA60_DECODED', u'DOY', u'QA_HOLLST',
            u'QA_FSEV1', u'NBR1', u'NBR2', u'STI', u'NDTI', u'CRC', u'REIP', u'GCVI', u'RDGCVI1', u'RDGCVI2',
            u'MTCI', u'MTCI2', u'WDRVI', u'GRWDRVI', u'RDWDRVI', u'RDNDVI1', u'RDNDVI2', u'NDVI'
        ]
        compare_bands(self, testimg_s2qa, expected_bands, {'msg': 'Sentinel 2 with VIs and cloud masks had the wrong bands'})


if __name__ == '__main__':
    unittest.main()
