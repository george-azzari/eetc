"""
# TODO: look into how to make this a unittest proper. Does CI work with EE?
"""
import unittest

import ee

from gee_tools.seasonal import getS1Plus, getS2Plus, sentinel_combined_medians
from tests.test_utils import compare_bands


class SeasonalTestCase(unittest.TestCase):

    def setUp(self):
        ee.Initialize()

        # ROI around Nairobi, Kenya.
        self.roi = ee.Geometry.Rectangle([
            36.7697414344982, -1.3481849306508715,
            37.095211405201326, -1.1278233090907424
        ])

    def test_s1(self):
        s1seascoll = getS1Plus(self.roi, 2018, True, True)
        s1seasimg = ee.Image(s1seascoll.first())

        expected_bands = [
            u'VV', u'VH', u'angle', u'LIA', u'VH_RLSPCK',
            u'VV_RLSPCK', u'DIFF', u'RATIO', u'DIFF_RLSPCK',
            u'RATIO_RLSPCK', u'DOY'
        ]
        compare_bands(self, s1seasimg, expected_bands, {'msg': 'Test image for S1 had the wrong bands'})

    def test_s2(self):
        s2seascoll = getS2Plus(self.roi, 2018, True)
        s2seasimg = ee.Image(s2seascoll.first())

        expected_bands = [
            u'AEROS', u'BLUE', u'GREEN', u'RED', u'RDED1', u'RDED2', u'RDED3', u'NIR', u'RDED4', u'VAPOR',
            u'CIRRU', u'SWIR1', u'SWIR2', u'QA10', u'QA20', u'QA60', u'QA60_DECODED', u'DOY', u'QA_HOLLST',
            u'QA_FSEV1', u'NBR1', u'NBR2', u'STI', u'NDTI', u'CRC', u'REIP', u'GCVI', u'RDGCVI1', u'RDGCVI2',
            u'MTCI', u'MTCI2', u'WDRVI', u'GRWDRVI', u'RDWDRVI', u'RDNDVI1', u'RDNDVI2', u'NDVI'
        ]
        compare_bands(self, s2seasimg, expected_bands, {'msg': 'Test image for S2 had the wrong bands'})

    def test_seasonal(self):
        glcmvars = ['contrast', 'corr', 'var', 'savg', 'prom']
        seasmed = sentinel_combined_medians(self.roi, 2018, False, False, True, True, glcmvars)
        # TODO can we generate this in code.
        expected_bands = [
            u'AEROS_S1', u'BLUE_S1', u'GREEN_S1', u'RED_S1', u'RDED1_S1', u'RDED2_S1', u'RDED3_S1', u'NIR_S1', u'RDED4_S1',
            u'VAPOR_S1', u'CIRRU_S1', u'SWIR1_S1', u'SWIR2_S1', u'QA10_S1', u'QA20_S1', u'QA60_S1', u'QA60_DECODED_S1',
            u'DOY_S1', u'QA_HOLLST_S1', u'QA_FSEV1_S1', u'AEROS_S2', u'BLUE_S2', u'GREEN_S2', u'RED_S2', u'RDED1_S2',
            u'RDED2_S2', u'RDED3_S2', u'NIR_S2', u'RDED4_S2', u'VAPOR_S2', u'CIRRU_S2', u'SWIR1_S2', u'SWIR2_S2', u'QA10_S2',
            u'QA20_S2', u'QA60_S2', u'QA60_DECODED_S2', u'DOY_S2', u'QA_HOLLST_S2', u'QA_FSEV1_S2', u'AEROS_S3', u'BLUE_S3',
            u'GREEN_S3', u'RED_S3', u'RDED1_S3', u'RDED2_S3', u'RDED3_S3', u'NIR_S3', u'RDED4_S3', u'VAPOR_S3', u'CIRRU_S3',
            u'SWIR1_S3', u'SWIR2_S3', u'QA10_S3', u'QA20_S3', u'QA60_S3', u'QA60_DECODED_S3', u'DOY_S3', u'QA_HOLLST_S3',
            u'QA_FSEV1_S3', u'AEROS_S1_contrast', u'BLUE_S1_contrast', u'GREEN_S1_contrast', u'RED_S1_contrast', u'RDED1_S1_contrast',
            u'RDED2_S1_contrast', u'RDED3_S1_contrast', u'NIR_S1_contrast', u'RDED4_S1_contrast', u'VAPOR_S1_contrast',
            u'CIRRU_S1_contrast', u'SWIR1_S1_contrast', u'SWIR2_S1_contrast', u'QA10_S1_contrast', u'QA20_S1_contrast',
            u'QA60_S1_contrast', u'QA60_DECODED_S1_contrast', u'DOY_S1_contrast', u'QA_HOLLST_S1_contrast', u'QA_FSEV1_S1_contrast',
            u'AEROS_S2_contrast', u'BLUE_S2_contrast', u'GREEN_S2_contrast', u'RED_S2_contrast', u'RDED1_S2_contrast', u'RDED2_S2_contrast',
            u'RDED3_S2_contrast', u'NIR_S2_contrast', u'RDED4_S2_contrast', u'VAPOR_S2_contrast', u'CIRRU_S2_contrast', u'SWIR1_S2_contrast',
            u'SWIR2_S2_contrast', u'QA10_S2_contrast', u'QA20_S2_contrast', u'QA60_S2_contrast', u'QA60_DECODED_S2_contrast',
            u'DOY_S2_contrast', u'QA_HOLLST_S2_contrast', u'QA_FSEV1_S2_contrast', u'AEROS_S3_contrast', u'BLUE_S3_contrast',
            u'GREEN_S3_contrast', u'RED_S3_contrast', u'RDED1_S3_contrast', u'RDED2_S3_contrast', u'RDED3_S3_contrast', u'NIR_S3_contrast',
            u'RDED4_S3_contrast', u'VAPOR_S3_contrast', u'CIRRU_S3_contrast', u'SWIR1_S3_contrast', u'SWIR2_S3_contrast', u'QA10_S3_contrast',
            u'QA20_S3_contrast', u'QA60_S3_contrast', u'QA60_DECODED_S3_contrast', u'DOY_S3_contrast', u'QA_HOLLST_S3_contrast', u'QA_FSEV1_S3_contrast',
            u'AEROS_S1_corr', u'BLUE_S1_corr', u'GREEN_S1_corr', u'RED_S1_corr', u'RDED1_S1_corr', u'RDED2_S1_corr', u'RDED3_S1_corr', u'NIR_S1_corr',
            u'RDED4_S1_corr', u'VAPOR_S1_corr', u'CIRRU_S1_corr', u'SWIR1_S1_corr', u'SWIR2_S1_corr', u'QA10_S1_corr', u'QA20_S1_corr',
            u'QA60_S1_corr', u'QA60_DECODED_S1_corr', u'DOY_S1_corr', u'QA_HOLLST_S1_corr', u'QA_FSEV1_S1_corr', u'AEROS_S2_corr',
            u'BLUE_S2_corr', u'GREEN_S2_corr', u'RED_S2_corr', u'RDED1_S2_corr', u'RDED2_S2_corr', u'RDED3_S2_corr', u'NIR_S2_corr',
            u'RDED4_S2_corr', u'VAPOR_S2_corr', u'CIRRU_S2_corr', u'SWIR1_S2_corr', u'SWIR2_S2_corr', u'QA10_S2_corr', u'QA20_S2_corr',
            u'QA60_S2_corr', u'QA60_DECODED_S2_corr', u'DOY_S2_corr', u'QA_HOLLST_S2_corr', u'QA_FSEV1_S2_corr', u'AEROS_S3_corr', u'BLUE_S3_corr',
            u'GREEN_S3_corr', u'RED_S3_corr', u'RDED1_S3_corr', u'RDED2_S3_corr', u'RDED3_S3_corr', u'NIR_S3_corr', u'RDED4_S3_corr', u'VAPOR_S3_corr',
            u'CIRRU_S3_corr', u'SWIR1_S3_corr', u'SWIR2_S3_corr', u'QA10_S3_corr', u'QA20_S3_corr', u'QA60_S3_corr', u'QA60_DECODED_S3_corr',
            u'DOY_S3_corr', u'QA_HOLLST_S3_corr', u'QA_FSEV1_S3_corr', u'AEROS_S1_var', u'BLUE_S1_var', u'GREEN_S1_var', u'RED_S1_var',
            u'RDED1_S1_var', u'RDED2_S1_var', u'RDED3_S1_var', u'NIR_S1_var', u'RDED4_S1_var', u'VAPOR_S1_var', u'CIRRU_S1_var', u'SWIR1_S1_var',
            u'SWIR2_S1_var', u'QA10_S1_var', u'QA20_S1_var', u'QA60_S1_var', u'QA60_DECODED_S1_var', u'DOY_S1_var', u'QA_HOLLST_S1_var', u'QA_FSEV1_S1_var',
            u'AEROS_S2_var', u'BLUE_S2_var', u'GREEN_S2_var', u'RED_S2_var', u'RDED1_S2_var', u'RDED2_S2_var', u'RDED3_S2_var', u'NIR_S2_var',
            u'RDED4_S2_var', u'VAPOR_S2_var', u'CIRRU_S2_var', u'SWIR1_S2_var', u'SWIR2_S2_var', u'QA10_S2_var', u'QA20_S2_var', u'QA60_S2_var',
            u'QA60_DECODED_S2_var', u'DOY_S2_var', u'QA_HOLLST_S2_var', u'QA_FSEV1_S2_var', u'AEROS_S3_var', u'BLUE_S3_var', u'GREEN_S3_var',
            u'RED_S3_var', u'RDED1_S3_var', u'RDED2_S3_var', u'RDED3_S3_var', u'NIR_S3_var', u'RDED4_S3_var', u'VAPOR_S3_var', u'CIRRU_S3_var',
            u'SWIR1_S3_var', u'SWIR2_S3_var', u'QA10_S3_var', u'QA20_S3_var', u'QA60_S3_var', u'QA60_DECODED_S3_var', u'DOY_S3_var', u'QA_HOLLST_S3_var',
            u'QA_FSEV1_S3_var', u'AEROS_S1_savg', u'BLUE_S1_savg', u'GREEN_S1_savg', u'RED_S1_savg', u'RDED1_S1_savg', u'RDED2_S1_savg', u'RDED3_S1_savg',
            u'NIR_S1_savg', u'RDED4_S1_savg', u'VAPOR_S1_savg', u'CIRRU_S1_savg', u'SWIR1_S1_savg', u'SWIR2_S1_savg', u'QA10_S1_savg', u'QA20_S1_savg',
            u'QA60_S1_savg', u'QA60_DECODED_S1_savg', u'DOY_S1_savg', u'QA_HOLLST_S1_savg', u'QA_FSEV1_S1_savg', u'AEROS_S2_savg', u'BLUE_S2_savg',
            u'GREEN_S2_savg', u'RED_S2_savg', u'RDED1_S2_savg', u'RDED2_S2_savg', u'RDED3_S2_savg', u'NIR_S2_savg', u'RDED4_S2_savg', u'VAPOR_S2_savg',
            u'CIRRU_S2_savg', u'SWIR1_S2_savg', u'SWIR2_S2_savg', u'QA10_S2_savg', u'QA20_S2_savg', u'QA60_S2_savg', u'QA60_DECODED_S2_savg', u'DOY_S2_savg',
            u'QA_HOLLST_S2_savg', u'QA_FSEV1_S2_savg', u'AEROS_S3_savg', u'BLUE_S3_savg', u'GREEN_S3_savg', u'RED_S3_savg', u'RDED1_S3_savg',
            u'RDED2_S3_savg', u'RDED3_S3_savg', u'NIR_S3_savg', u'RDED4_S3_savg', u'VAPOR_S3_savg', u'CIRRU_S3_savg', u'SWIR1_S3_savg', u'SWIR2_S3_savg',
            u'QA10_S3_savg', u'QA20_S3_savg', u'QA60_S3_savg', u'QA60_DECODED_S3_savg', u'DOY_S3_savg', u'QA_HOLLST_S3_savg', u'QA_FSEV1_S3_savg',
            u'AEROS_S1_prom', u'BLUE_S1_prom', u'GREEN_S1_prom', u'RED_S1_prom', u'RDED1_S1_prom', u'RDED2_S1_prom', u'RDED3_S1_prom',
            u'NIR_S1_prom', u'RDED4_S1_prom', u'VAPOR_S1_prom', u'CIRRU_S1_prom', u'SWIR1_S1_prom', u'SWIR2_S1_prom', u'QA10_S1_prom',
            u'QA20_S1_prom', u'QA60_S1_prom', u'QA60_DECODED_S1_prom', u'DOY_S1_prom', u'QA_HOLLST_S1_prom', u'QA_FSEV1_S1_prom', u'AEROS_S2_prom',
            u'BLUE_S2_prom', u'GREEN_S2_prom', u'RED_S2_prom', u'RDED1_S2_prom', u'RDED2_S2_prom', u'RDED3_S2_prom', u'NIR_S2_prom',
            u'RDED4_S2_prom', u'VAPOR_S2_prom', u'CIRRU_S2_prom', u'SWIR1_S2_prom', u'SWIR2_S2_prom', u'QA10_S2_prom', u'QA20_S2_prom',
            u'QA60_S2_prom', u'QA60_DECODED_S2_prom', u'DOY_S2_prom', u'QA_HOLLST_S2_prom', u'QA_FSEV1_S2_prom', u'AEROS_S3_prom', u'BLUE_S3_prom',
            u'GREEN_S3_prom', u'RED_S3_prom', u'RDED1_S3_prom', u'RDED2_S3_prom', u'RDED3_S3_prom', u'NIR_S3_prom', u'RDED4_S3_prom', u'VAPOR_S3_prom',
            u'CIRRU_S3_prom', u'SWIR1_S3_prom', u'SWIR2_S3_prom', u'QA10_S3_prom', u'QA20_S3_prom', u'QA60_S3_prom', u'QA60_DECODED_S3_prom', u'DOY_S3_prom',
            u'QA_HOLLST_S3_prom', u'QA_FSEV1_S3_prom', u'VV_S1', u'VH_S1', u'angle_S1', u'VH_RLSPCK_S1', u'VV_RLSPCK_S1', u'DIFF_S1', u'RATIO_S1',
            u'DIFF_RLSPCK_S1', u'RATIO_RLSPCK_S1', u'DOY_S1_1', u'VV_S2', u'VH_S2', u'angle_S2', u'VH_RLSPCK_S2', u'VV_RLSPCK_S2', u'DIFF_S2', u'RATIO_S2',
            u'DIFF_RLSPCK_S2', u'RATIO_RLSPCK_S2', u'DOY_S2_1', u'VV_S3', u'VH_S3', u'angle_S3', u'VH_RLSPCK_S3', u'VV_RLSPCK_S3', u'DIFF_S3', u'RATIO_S3',
            u'DIFF_RLSPCK_S3', u'RATIO_RLSPCK_S3', u'DOY_S3_1', u'VV_S1_contrast', u'VH_S1_contrast', u'angle_S1_contrast', u'VH_RLSPCK_S1_contrast',
            u'VV_RLSPCK_S1_contrast', u'DIFF_S1_contrast', u'RATIO_S1_contrast', u'DIFF_RLSPCK_S1_contrast', u'RATIO_RLSPCK_S1_contrast', u'DOY_S1_contrast_1',
            u'VV_S2_contrast', u'VH_S2_contrast', u'angle_S2_contrast', u'VH_RLSPCK_S2_contrast', u'VV_RLSPCK_S2_contrast',
            u'DIFF_S2_contrast', u'RATIO_S2_contrast',
            u'DIFF_RLSPCK_S2_contrast', u'RATIO_RLSPCK_S2_contrast', u'DOY_S2_contrast_1', u'VV_S3_contrast', u'VH_S3_contrast', u'angle_S3_contrast',
            u'VH_RLSPCK_S3_contrast', u'VV_RLSPCK_S3_contrast', u'DIFF_S3_contrast', u'RATIO_S3_contrast',
            u'DIFF_RLSPCK_S3_contrast', u'RATIO_RLSPCK_S3_contrast',
            u'DOY_S3_contrast_1', u'VV_S1_corr', u'VH_S1_corr', u'angle_S1_corr', u'VH_RLSPCK_S1_corr', u'VV_RLSPCK_S1_corr', u'DIFF_S1_corr', u'RATIO_S1_corr',
            u'DIFF_RLSPCK_S1_corr', u'RATIO_RLSPCK_S1_corr', u'DOY_S1_corr_1', u'VV_S2_corr', u'VH_S2_corr', u'angle_S2_corr', u'VH_RLSPCK_S2_corr',
            u'VV_RLSPCK_S2_corr', u'DIFF_S2_corr', u'RATIO_S2_corr', u'DIFF_RLSPCK_S2_corr', u'RATIO_RLSPCK_S2_corr', u'DOY_S2_corr_1', u'VV_S3_corr',
            u'VH_S3_corr', u'angle_S3_corr', u'VH_RLSPCK_S3_corr', u'VV_RLSPCK_S3_corr', u'DIFF_S3_corr', u'RATIO_S3_corr', u'DIFF_RLSPCK_S3_corr',
            u'RATIO_RLSPCK_S3_corr', u'DOY_S3_corr_1', u'VV_S1_var', u'VH_S1_var', u'angle_S1_var', u'VH_RLSPCK_S1_var', u'VV_RLSPCK_S1_var', u'DIFF_S1_var',
            u'RATIO_S1_var', u'DIFF_RLSPCK_S1_var', u'RATIO_RLSPCK_S1_var', u'DOY_S1_var_1', u'VV_S2_var', u'VH_S2_var', u'angle_S2_var', u'VH_RLSPCK_S2_var',
            u'VV_RLSPCK_S2_var', u'DIFF_S2_var', u'RATIO_S2_var', u'DIFF_RLSPCK_S2_var', u'RATIO_RLSPCK_S2_var', u'DOY_S2_var_1', u'VV_S3_var', u'VH_S3_var',
            u'angle_S3_var', u'VH_RLSPCK_S3_var', u'VV_RLSPCK_S3_var', u'DIFF_S3_var', u'RATIO_S3_var', u'DIFF_RLSPCK_S3_var', u'RATIO_RLSPCK_S3_var',
            u'DOY_S3_var_1', u'VV_S1_savg', u'VH_S1_savg', u'angle_S1_savg', u'VH_RLSPCK_S1_savg', u'VV_RLSPCK_S1_savg', u'DIFF_S1_savg', u'RATIO_S1_savg',
            u'DIFF_RLSPCK_S1_savg', u'RATIO_RLSPCK_S1_savg', u'DOY_S1_savg_1', u'VV_S2_savg', u'VH_S2_savg', u'angle_S2_savg', u'VH_RLSPCK_S2_savg',
            u'VV_RLSPCK_S2_savg', u'DIFF_S2_savg', u'RATIO_S2_savg', u'DIFF_RLSPCK_S2_savg', u'RATIO_RLSPCK_S2_savg', u'DOY_S2_savg_1', u'VV_S3_savg',
            u'VH_S3_savg', u'angle_S3_savg', u'VH_RLSPCK_S3_savg', u'VV_RLSPCK_S3_savg', u'DIFF_S3_savg', u'RATIO_S3_savg', u'DIFF_RLSPCK_S3_savg',
            u'RATIO_RLSPCK_S3_savg', u'DOY_S3_savg_1', u'VV_S1_prom', u'VH_S1_prom', u'angle_S1_prom', u'VH_RLSPCK_S1_prom',
            u'VV_RLSPCK_S1_prom', u'DIFF_S1_prom',
            u'RATIO_S1_prom', u'DIFF_RLSPCK_S1_prom', u'RATIO_RLSPCK_S1_prom', u'DOY_S1_prom_1', u'VV_S2_prom', u'VH_S2_prom', u'angle_S2_prom',
            u'VH_RLSPCK_S2_prom', u'VV_RLSPCK_S2_prom', u'DIFF_S2_prom', u'RATIO_S2_prom', u'DIFF_RLSPCK_S2_prom', u'RATIO_RLSPCK_S2_prom', u'DOY_S2_prom_1',
            u'VV_S3_prom', u'VH_S3_prom', u'angle_S3_prom', u'VH_RLSPCK_S3_prom', u'VV_RLSPCK_S3_prom', u'DIFF_S3_prom',
            u'RATIO_S3_prom', u'DIFF_RLSPCK_S3_prom', u'RATIO_RLSPCK_S3_prom', u'DOY_S3_prom_1'
        ]
        compare_bands(self, seasmed, expected_bands, {'msg': 'Test image for sentinel_combined_medians had the wrong bands'})


if __name__ == '__main__':
    unittest.main()
