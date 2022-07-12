"""
python -m tests.datasources.test_optical_datasources
"""
import unittest

import ee

from gee_tools.datasources import optical_datasources as opt_ds


class LandsatSRTestCase(unittest.TestCase):

    def setUp(self):
        ee.Initialize()

        # ROI around Nairobi, Kenya.
        self.roi = ee.Geometry.Rectangle([
            36.7697414344982, -1.3481849306508715,
            37.095211405201326, -1.1278233090907424
        ])

    def test_get_quality_pixel_count(self):
        quality_pixel_count = opt_ds.LandsatSR(
            self.roi, '2012-08-01', '2012-12-31',
        ).get_quality_pixel_count()

        quality_pixel_count_band_names = ee.Dictionary({
            k: v.bandNames() for k, v
            in quality_pixel_count.items()
        })
        quality_pixel_count_band_names = quality_pixel_count_band_names.getInfo()

        self.assertDictEqual(
            {k: ['count'] for k in quality_pixel_count_band_names.keys()},
            quality_pixel_count_band_names
        )

    def test_get_quality_pixel_count(self):
        quality_pixel_count = opt_ds.LandsatSRQuality(
            self.roi, '2012-08-01', '2012-12-31',
        ).get_img_coll().mosaic()

        quality_pixel_count = quality_pixel_count.getInfo()
        bands = quality_pixel_count['bands']
        bands = {band_info['id']: band_info['data_type']['precision'] for band_info in bands}

        self.assertDictEqual(
            {
                band: 'int'
                for band in ['LS5_QUAL', 'LS7_QUAL', 'LS8_QUAL', 'LS5_TOTAL', 'LS7_TOTAL', 'LS8_TOTAL']
            },
            bands
        )


if __name__ == '__main__':
    unittest.main()
