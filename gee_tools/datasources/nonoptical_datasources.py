"""
Author: Anthony Perez

Each class represents a datasource and will have a property containing an image collection.
"""
import ee
from gee_tools.datasources.interface import MultiImageDatasource, GlobalImageDatasource, SingleImageDatasource, DatasourceError


class NightlightDatasource(GlobalImageDatasource):
    """Abstract class for nightlights datasources"""

    def init_coll(self, name):
        return ee.ImageCollection(name) \
                 .filterDate(self.start_date, self.end_date) \
                 .map(self.nl_rename) \
                 .sort('system:time_start')

    @staticmethod
    def nl_rename(scene):
        """scene.select([0], ['NIGHTLIGHTS'])"""
        return scene.select([0], ['NIGHTLIGHTS'])


class DMSPUncal(NightlightDatasource):
    """
    Uncalibrated DMSP nightlights
    Data in property dmsp
    """

    def build_img_coll(self):
        self.dmsp = self.init_coll("NOAA/DMSP-OLS/NIGHTTIME_LIGHTS")

    def get_img_coll(self):
        return self.dmsp


class DMSPCalV4(NightlightDatasource):
    """
    Calibrated DMSP nightlights
    Data in property dmsp
    """

    def build_img_coll(self):
        self.dmsp = self.init_coll("NOAA/DMSP-OLS/CALIBRATED_LIGHTS_V4")

    def get_img_coll(self):
        return self.dmsp


class VIIRSMonthlyStrCorr(NightlightDatasource):
    """
    Calibrated VIIRS nightlights
    Data in property viirs
    """

    def build_img_coll(self):
        self.viirs = self.init_coll("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG")

    def get_img_coll(self):
        return self.viirs


class DMSPCalVIIRSJoined(NightlightDatasource):
    """
    Returns the VIIRS image collection 'NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG'
    if the request date range is either after 2014 or contains 2014.
    2014 is the starting year for the VIIRS image collection listed above.

    Otherwise returns the calibrated DMSP image collection "NOAA/DMSP-OLS/CALIBRATED_LIGHTS_V4"
    """

    def build_img_coll(self):
        viirs = self.init_coll("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG")
        dmsp = self.init_coll("NOAA/DMSP-OLS/CALIBRATED_LIGHTS_V4")

        # Carefull not to communicate to the server
        requested_daterange = ee.DateRange(self.start_date, self.end_date)
        viirs_start = ee.Date('2014-1-1')
        use_viirs = requested_daterange.contains(viirs_start)

        start_year = ee.Number(ee.Date(self.start_date).get('year'))
        viirs_start_year = ee.Number(viirs_start.get('year'))
        use_viirs2 = start_year.gte(viirs_start_year)

        self.nl = ee.Algorithms.If(use_viirs, viirs, dmsp)
        self.nl = ee.Algorithms.If(use_viirs2, viirs, self.nl)
        self.nl = ee.ImageCollection(self.nl)
    
    def get_img_coll(self):
        # TODO: not sure why this is helpful. All it's doing is returning a class method (George).
        return self.nl


class SRTMElevation(SingleImageDatasource):
    """
    Note: Near global
    """

    def build_img_coll(self):
        topo = ee.Image("USGS/SRTMGL1_003")
        band_names = ['ELEV', 'SLO', 'ASP']
        self.topo = ee.Algorithms.Terrain(topo).select(['elevation', 'slope', 'aspect'], band_names)
        self.coll = ee.ImageCollection([self.topo])

    def get_img_coll(self):
        return self.coll


class Palsar(GlobalImageDatasource):
    """
    Image Properties

    HH	    HH polarization backscattering coefficient, 16-bit DN.
    HV	    HV polarization backscattering coefficient, 16-bit DN.
    angle	Local incidence angle (degrees).
    date	Observation date (days since Jan 1, 1970).
    qa	    Processing information.

    Renamed to ["HH", "HV", "ANGLE", "DATE", "QA"]
    """

    def build_img_coll(self):

        self.palsar = ee.ImageCollection("JAXA/ALOS/PALSAR/YEARLY/SAR") \
                        .filterDate(self.start_date, self.end_date) \
                        .map(self.rename_pulsar) \
                        .map(self.mask_qa) \
                        .sort('system:time_start')

    def get_img_coll(self):
        return self.palsar

    @staticmethod
    def rename_pulsar(scene):
        """
        Image Properties

        HH	    HH polarization backscattering coefficient, 16-bit DN.
        HV	    HV polarization backscattering coefficient, 16-bit DN.
        angle	Local incidence angle (degrees).
        date	Observation date (days since Jan 1, 1970).
        qa	    Processing information.
        """
        band_names = ["HH", "HV", "ANGLE", "DATE", "QA"]
        return scene.select(range(len(band_names)), band_names)

    @staticmethod
    def decode_qa(scene):
        """
        Value	Color	Description
        0	    000000	No data
        50	    0000FF	Ocean and water
        100	    AAAA00	Radar layover
        150	    005555	Radar shadowing
        255	    AA9988	Land
        """
        qa = scene.select(["QA"])

        nodata = qa.eq(0)
        nodata = nodata.updateMask(nodata).rename(["pxqa_nodata"])

        hasdata = qa.neq(0)
        hasdata = hasdata.updateMask(hasdata).rename(["pxqa_hasdata"])

        water = qa.eq(50)
        water = water.updateMask(water).rename(["pxqa_water"])

        radar_layover = qa.eq(100)
        radar_layover = radar_layover.updateMask(radar_layover).rename(["pxqa_radar_layover"])

        radar_shad = qa.eq(150)
        radar_shad = radar_shad.updateMask(radar_shad).rename(["pxqa_radar_shad"])

        land = qa.eq(255)
        land = land.updateMask(land).rename(["pxqa_land"])

        masks = ee.Image.cat([nodata, hasdata, water, radar_layover, radar_shad, land])
        return masks

    def mask_qa(self, scene):
        masks = self.decode_qa(scene)
        scene.updateMask(masks.select(["pxqa_hasdata"]))
        return scene
