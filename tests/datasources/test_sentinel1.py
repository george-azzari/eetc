"""
# TODO: look into how to make this a unittest proper. Does CI work with EE?
"""
import unittest

import ee

from gee_tools.datasources import sentinel1
from tests.test_utils import compare_bands


class Sentinel1TestCase(unittest.TestCase):

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
        Test base collection (raw S1 assets).
        """
        # Generate base collection.
        s1base = sentinel1.Sentinel1(self.roi, self.sdate, self.edate)
        testimg_base = ee.Image(s1base.coll.first())  # TODO tests should be agnostic to the underlying implementation

        expected_bands = [u'VV', u'VH', u'angle']
        compare_bands(self, testimg_base, expected_bands, {'msg': 'Sentinel 1 base image had the wrong bands'})

    def test_extra_bands(self):
        """
        Test extra bands (difference and polarization).
        """
        # Get collection with extra bands.
        s1base = sentinel1.Sentinel1(self.roi, self.sdate, self.edate)
        s1extra = s1base.get_img_coll(correctlia=False, addbands=True, addspeckle=False, addtexture=False, orbit='ascending')
        testimg_s1extra = ee.Image(s1extra.first())

        expected_bands = [u'VV', u'VH', u'angle', u'DIFF', u'RATIO']
        compare_bands(self, testimg_s1extra, expected_bands, {'msg': 'Sentinel 1 with extras had the wrong bands'})
    
    def test_extra_bands_constructor(self):
        """
        Test extra bands (difference and polarization).
        """
        # Get collection with extra bands.
        s1base = sentinel1.Sentinel1(
            self.roi, self.sdate, self.edate,
            correctlia=False, addbands=True,
            addspeckle=False, addtexture=False,
            orbit='ascending',
        )
        s1extra = s1base.get_img_coll()
        testimg_s1extra = ee.Image(s1extra.first())

        expected_bands = [u'VV', u'VH', u'angle', u'DIFF', u'RATIO']
        compare_bands(self, testimg_s1extra, expected_bands, {'msg': 'Sentinel 1 with extras had the wrong bands'})

    def test_speckle_correction(self):
        """
        Test speckle-correction.
        """
        # Get collection with extra bands.
        s1base = sentinel1.Sentinel1(self.roi, self.sdate, self.edate)
        s1speckle = s1base.get_img_coll(correctlia=False, addbands=True, addspeckle=True, addtexture=False, orbit='ascending')
        testimg_s1speckle = ee.Image(s1speckle.first())

        expected_bands = [u'VV', u'VH', u'angle', u'VH_RLSPCK', u'VV_RLSPCK', u'DIFF', u'RATIO', u'DIFF_RLSPCK', u'RATIO_RLSPCK']
        compare_bands(self, testimg_s1speckle, expected_bands, {'msg': 'Sentinel 1 with speckle-correction had the wrong bands'})

    def test_speckle_correction_constructor(self):
        """
        Test speckle-correction.
        """
        # Get collection with extra bands.
        s1base = sentinel1.Sentinel1(
            self.roi, self.sdate, self.edate,
            correctlia=False, addbands=True,
            addspeckle=True, addtexture=False,
            orbit='ascending',
        )
        s1speckle = s1base.get_img_coll()
        testimg_s1speckle = ee.Image(s1speckle.first())

        expected_bands = [u'VV', u'VH', u'angle', u'VH_RLSPCK', u'VV_RLSPCK', u'DIFF', u'RATIO', u'DIFF_RLSPCK', u'RATIO_RLSPCK']
        compare_bands(self, testimg_s1speckle, expected_bands, {'msg': 'Sentinel 1 with speckle-correction had the wrong bands'})

    def test_glcm_texture(self):
        """
        Test GLCM texture.
        """
        # Get collection with extra bands.
        s1base = sentinel1.Sentinel1(self.roi, self.sdate, self.edate)
        s1glcm = s1base.get_img_coll(correctlia=False, addbands=True, addspeckle=True, addtexture=True, orbit='ascending')
        testimg_s1glcm = ee.Image(s1glcm.first())

        # TODO can we generate this with code
        # TODO not sure if this is correct
        expected_bands = [
            u'VV', u'VH', u'angle', u'VH_RLSPCK', u'VV_RLSPCK', u'DIFF', u'RATIO', u'DIFF_RLSPCK', u'RATIO_RLSPCK', u'DIFF_asm',
            u'DIFF_contrast', u'DIFF_corr', u'DIFF_var', u'DIFF_idm', u'DIFF_savg', u'DIFF_svar', u'DIFF_sent', u'DIFF_ent',
            u'DIFF_dvar', u'DIFF_dent', u'DIFF_imcorr1', u'DIFF_imcorr2', u'DIFF_maxcorr', u'DIFF_diss', u'DIFF_inertia',
            u'DIFF_shade', u'DIFF_prom', u'DIFF_RLSPCK_asm', u'DIFF_RLSPCK_contrast', u'DIFF_RLSPCK_corr', u'DIFF_RLSPCK_var',
            u'DIFF_RLSPCK_idm', u'DIFF_RLSPCK_savg', u'DIFF_RLSPCK_svar', u'DIFF_RLSPCK_sent', u'DIFF_RLSPCK_ent', u'DIFF_RLSPCK_dvar',
            u'DIFF_RLSPCK_dent', u'DIFF_RLSPCK_imcorr1', u'DIFF_RLSPCK_imcorr2', u'DIFF_RLSPCK_maxcorr', u'DIFF_RLSPCK_diss',
            u'DIFF_RLSPCK_inertia', u'DIFF_RLSPCK_shade', u'DIFF_RLSPCK_prom', u'VH_asm', u'VH_contrast', u'VH_corr', u'VH_var',
            u'VH_idm', u'VH_savg', u'VH_svar', u'VH_sent', u'VH_ent', u'VH_dvar', u'VH_dent', u'VH_imcorr1', u'VH_imcorr2', u'VH_maxcorr',
            u'VH_diss', u'VH_inertia', u'VH_shade', u'VH_prom', u'VH_RLSPCK_asm', u'VH_RLSPCK_contrast', u'VH_RLSPCK_corr',
            u'VH_RLSPCK_var', u'VH_RLSPCK_idm', u'VH_RLSPCK_savg', u'VH_RLSPCK_svar', u'VH_RLSPCK_sent', u'VH_RLSPCK_ent',
            u'VH_RLSPCK_dvar', u'VH_RLSPCK_dent', u'VH_RLSPCK_imcorr1', u'VH_RLSPCK_imcorr2', u'VH_RLSPCK_maxcorr', u'VH_RLSPCK_diss',
            u'VH_RLSPCK_inertia', u'VH_RLSPCK_shade', u'VH_RLSPCK_prom', u'VV_asm', u'VV_contrast', u'VV_corr', u'VV_var', u'VV_idm',
            u'VV_savg', u'VV_svar', u'VV_sent', u'VV_ent', u'VV_dvar', u'VV_dent', u'VV_imcorr1', u'VV_imcorr2', u'VV_maxcorr',
            u'VV_diss', u'VV_inertia', u'VV_shade', u'VV_prom', u'VV_RLSPCK_asm', u'VV_RLSPCK_contrast', u'VV_RLSPCK_corr',
            u'VV_RLSPCK_var', u'VV_RLSPCK_idm', u'VV_RLSPCK_savg', u'VV_RLSPCK_svar', u'VV_RLSPCK_sent', u'VV_RLSPCK_ent',
            u'VV_RLSPCK_dvar', u'VV_RLSPCK_dent', u'VV_RLSPCK_imcorr1', u'VV_RLSPCK_imcorr2', u'VV_RLSPCK_maxcorr', u'VV_RLSPCK_diss',
            u'VV_RLSPCK_inertia', u'VV_RLSPCK_shade', u'VV_RLSPCK_prom'
        ]
        compare_bands(self, testimg_s1glcm, expected_bands, {'msg': 'Sentinel 1 with speckle-correction and GLCM texture had the wrong bands'})


if __name__ == '__main__':
    unittest.main()
