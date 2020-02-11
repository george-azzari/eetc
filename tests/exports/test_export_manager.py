"""
Author: Anthony Perez

python -m tests.exports.test_export_manager
"""
import unittest
import itertools

import ee

from gee_tools.datasources.optical_datasources import LandsatSR, MODISnbar
from gee_tools.datasources.pop_and_urban_datasources import GHSLPop, CityAccessibility
from gee_tools.exports.export_manager import ExportManager
from gee_tools.exports.constants import EPSG3857


IN_MODIS_BANDS = [u'RED', u'NIR', u'BLUE', u'GREEN', u'SWIR1', u'SWIR2']
OUT_MODIS_BANDS = [u'MODIS_{}'.format(b) for b in IN_MODIS_BANDS]
IN_LX_BANDS = [u'BLUE', u'GREEN', u'RED', u'NIR', u'SWIR1', u'SWIR2', u'TEMP1']
OUT_LX_BANDS = [u'LX_{}'.format(b) for b in IN_LX_BANDS]


def modis_composite(img_coll):
    return img_coll.select(
        IN_MODIS_BANDS, OUT_MODIS_BANDS
    ).median()


def lx_composite(img_coll):
    return img_coll.select(
        IN_LX_BANDS, OUT_LX_BANDS
    ).median()


class ExportManagerTestCase(unittest.TestCase):

    def setUp(self):
        ee.Initialize()

    def test_base_collection(self):
        """
        Test basic configuration.
        """
        config = {
            "landsat": {
                "class": LandsatSR,
                "args": {},
                "composite_fn": lx_composite,
                "bands": OUT_LX_BANDS,
            },
            'modis': {
                "class": MODISnbar,
                "args": {},
                "composite_fn": modis_composite,
                "bands": OUT_MODIS_BANDS,
            },
            'ghsl_pop': {
                "class": GHSLPop,
                "args": { 'use_closest_image': True },
                "composite_fn": lambda img_coll: img_coll.mean(),
                "bands": ['POPULATION'],
            },
        }
        config_with_lat = {
            "landsat": {
                "class": LandsatSR,
                "args": {},
                "composite_fn": lx_composite,
                "bands": OUT_LX_BANDS,
            },
            'modis': {
                "class": MODISnbar,
                "args": {},
                "composite_fn": modis_composite,
                "bands": OUT_MODIS_BANDS,
            },
            'ghsl_pop': {
                "class": GHSLPop,
                "args": { 'use_closest_image': True },
                "composite_fn": lambda img_coll: img_coll.mean().addBands(
                    ee.Image(0).select([0], ['LAT'])
                ),
                "bands": ['POPULATION', 'LAT'],
            },
        }

        image_spec = {
            'start_date': '2011-01-01',
            'end_date': '2011-12-31',
            'filterpoly': ee.Geometry.Point(38.76, 9.01).buffer(3000).bounds(),
            'projection': EPSG3857,  # CRS
            'scale': 30,
        }

        export_manager = ExportManager(config)
        scene, scene_reported_output_bands = export_manager.get_scene(image_spec)

        test_fc = ee.FeatureCollection([
            ee.Feature(ee.Geometry.Point(38.76, 9.01), {})
        ])
        # Sample 255 by 255 tiles around the given points.
        sampled_fc, samples_reported_output_bands = export_manager.sample_tiles(test_fc, image_spec, 127)

        expected_bands = list(itertools.chain(
            *[c['bands'] for c in config.values()]
        )) + [u'LAT', u'LON']
        actual_bands = scene.bandNames().getInfo()

        self.assertEqual(sorted(expected_bands), sorted(actual_bands))
        self.assertEqual(sorted(expected_bands), sorted(scene_reported_output_bands))
        self.assertEqual(sorted(scene_reported_output_bands), sorted(actual_bands))
        self.assertEqual(sorted(scene_reported_output_bands), sorted(samples_reported_output_bands))

        sampled_fc = sampled_fc.getInfo()
        feat = sampled_fc['features'][0]
        props = feat['properties']

        for band in expected_bands:
            self.assertTrue(band in props)
            band_value = props[band]

            self.assertEqual(len(band_value), 255)
            self.assertEqual(len(band_value[0]), 255)

        export_manager = ExportManager(config_with_lat)
        with self.assertRaises(ValueError):
            export_manager.get_scene(image_spec)

    def test_empty_collection(self):
        """
        Test what happens when a collection is filtered to the empty collection.

        Empty collection are omitted.
        """
        config = {
            'modis': {
                "class": MODISnbar,
                "args": {},
                "composite_fn": modis_composite,
                "bands": OUT_MODIS_BANDS,
            },
            'accessibility': {
                "class": CityAccessibility,
                "args": {},
                "composite_fn": None,
                "bands": ['ACCESSIBILITY'],
            }
        }

        image_spec = {
            'start_date': '1950-01-01',
            'end_date': '1950-12-31',
            'filterpoly': ee.Geometry.Point(38.76, 9.01).buffer(3000).bounds(),
            'projection': EPSG3857,  # CRS
            'scale': 30,
        }

        export_manager = ExportManager(config)
        scene, scene_reported_output_bands = export_manager.get_scene(image_spec)

        test_fc = ee.FeatureCollection([
            ee.Feature(ee.Geometry.Point(38.76, 9.01), {})
        ])
        # Sample 255 by 255 tiles around the given points.
        sampled_fc, samples_reported_output_bands = export_manager.sample_tiles(test_fc, image_spec, 127)

        expected_bands = ['ACCESSIBILITY', 'LAT', 'LON']
        actual_bands = scene.bandNames().getInfo()

        self.assertEqual(sorted(expected_bands), sorted(actual_bands))
        self.assertEqual(
            sorted(samples_reported_output_bands),
            sorted(OUT_MODIS_BANDS + expected_bands)
        )

        sampled_fc = sampled_fc.getInfo()
        feat = sampled_fc['features'][0]
        props = feat['properties']

        for band in OUT_MODIS_BANDS:
            self.assertFalse(band in props)

    def test_empty_collection2(self):
        """
        Test what happens when a collection is filtered to the empty collection.

        Empty collection are omitted.
        """
        config = {
            'modis': {
                "class": MODISnbar,
                "args": {},
                "composite_fn": modis_composite,
                "bands": OUT_MODIS_BANDS,
            },
        }

        image_spec = {
            'start_date': '1950-01-01',
            'end_date': '1950-12-31',
            'filterpoly': ee.Geometry.Point(38.76, 9.01).buffer(3000).bounds(),
            'projection': EPSG3857,  # CRS
            'scale': 30,
        }

        export_manager = ExportManager(config)
        scene, scene_reported_output_bands = export_manager.get_scene(image_spec)

        test_fc = ee.FeatureCollection([
            ee.Feature(ee.Geometry.Point(38.76, 9.01), {})
        ])
        # Sample 255 by 255 tiles around the given points.
        sampled_fc, samples_reported_output_bands = export_manager.sample_tiles(test_fc, image_spec, 127)

        expected_bands = []
        actual_bands = scene.bandNames().getInfo()

        self.assertEqual(sorted(expected_bands), sorted(actual_bands))
        self.assertEqual(
            sorted(samples_reported_output_bands),
            sorted(OUT_MODIS_BANDS + ['LAT', 'LON'])
        )

        sampled_fc = sampled_fc.getInfo()
        feat = sampled_fc['features'][0]
        props = feat['properties']

        for band in OUT_MODIS_BANDS:
            self.assertFalse(band in props)

    def test_sample_tiles_unstacked(self):
        """
        Test basic configuration.
        """
        config = {
            "landsat": {
                "class": LandsatSR,
                "args": {},
                "composite_fn": lx_composite,
                "bands": OUT_LX_BANDS,
            },
            'modis': {
                "class": MODISnbar,
                "args": {},
                "composite_fn": modis_composite,
                "bands": OUT_MODIS_BANDS,
            },
            'ghsl_pop': {
                "class": GHSLPop,
                "args": { 'use_closest_image': True },
                "composite_fn": lambda img_coll: img_coll.mean(),
                "bands": ['POPULATION'],
            },
        }

        image_spec = {
            'start_date': '2011-01-01',
            'end_date': '2011-12-31',
            'filterpoly': ee.Geometry.Point(38.76, 9.01).buffer(3000).bounds(),
            'projection': EPSG3857,  # CRS
            'scale': 30,
        }

        export_manager = ExportManager(config)
        scene, scene_reported_output_bands = export_manager.get_scene(image_spec)

        test_fc = ee.FeatureCollection([
            ee.Feature(ee.Geometry.Point(38.76, 9.01), {})
        ])
        # Sample 255 by 255 tiles around the given points.
        sampled_fc, samples_reported_output_bands = export_manager.sample_tiles_unstacked(test_fc, image_spec, 127)

        expected_bands = list(itertools.chain(
            *[c['bands'] for c in config.values()]
        )) + [u'LAT', u'LON']
        actual_bands = scene.bandNames().getInfo()

        self.assertEqual(sorted(expected_bands), sorted(actual_bands))
        self.assertEqual(sorted(expected_bands), sorted(scene_reported_output_bands))
        self.assertEqual(sorted(scene_reported_output_bands), sorted(actual_bands))
        self.assertEqual(sorted(scene_reported_output_bands), sorted(samples_reported_output_bands))

        sampled_fc = sampled_fc.getInfo()
        feat = sampled_fc['features'][0]
        props = feat['properties']

        for band in expected_bands:
            self.assertTrue(band in props)
            band_value = props[band]

            self.assertEqual(len(band_value), 255)
            self.assertEqual(len(band_value[0]), 255)

        # Test included LAT band

        config_with_lat = {
            "landsat": {
                "class": LandsatSR,
                "args": {},
                "composite_fn": lx_composite,
                "bands": OUT_LX_BANDS,
            },
            'modis': {
                "class": MODISnbar,
                "args": {},
                "composite_fn": modis_composite,
                "bands": OUT_MODIS_BANDS,
            },
            'ghsl_pop': {
                "class": GHSLPop,
                "args": { 'use_closest_image': True },
                "composite_fn": lambda img_coll: img_coll.mean().addBands(
                    ee.Image(0).select([0], ['LAT'])
                ),
                "bands": ['POPULATION', 'LAT'],
            },
        }

        export_manager = ExportManager(config_with_lat)
        with self.assertRaises(ValueError):
            export_manager.sample_tiles_unstacked(test_fc, image_spec, 127)

        # Test empty collection

        config = {
            'modis': {
                "class": MODISnbar,
                "args": {},
                "composite_fn": modis_composite,
                "bands": OUT_MODIS_BANDS,
            },
        }

        image_spec = {
            'start_date': '1950-01-01',
            'end_date': '1950-12-31',
            'filterpoly': ee.Geometry.Point(38.76, 9.01).buffer(3000).bounds(),
            'projection': EPSG3857,  # CRS
            'scale': 30,
        }

        export_manager = ExportManager(config)
        sampled_fc, samples_reported_output_bands = export_manager.sample_tiles(test_fc, image_spec, 127)

        self.assertEqual(
            sorted(samples_reported_output_bands),
            sorted(OUT_MODIS_BANDS + ['LAT', 'LON'])
        )

        sampled_fc = sampled_fc.getInfo()
        feat = sampled_fc['features'][0]
        props = feat['properties']

        for band in OUT_MODIS_BANDS:
            self.assertFalse(band in props)


if __name__ == '__main__':
    unittest.main()
