"""
Author: Anthony Perez

Each class represents a datasource and will have a property containing an image collection.
"""
import datetime

import ee

from gee_tools.datasources.interface import MultiImageDatasource, GlobalImageDatasource, SingleImageDatasource, DatasourceError
from gee_tools.datasources.util import get_center_date, get_closest_to_date, compute_monthly_seasonal_annual_and_long_running


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


class DMSPCalV4ClosestBefore2014(GlobalImageDatasource):
    """
    The closest image in the 'NOAA/DMSP-OLS/CALIBRATED_LIGHTS_V4'
    collection to the given filter date.  If the given filter date includes
    January 1st 2014 or after, the output is all 0s.
    """
    @staticmethod
    def band_names():
        return ['DMSP']

    def build_img_coll(self):
        dmsp_original = ee.FeatureCollection("NOAA/DMSP-OLS/CALIBRATED_LIGHTS_V4")

        center_date = get_center_date(self.start_date, self.end_date)
        dmsp_img = get_closest_to_date(dmsp_original, center_date)
        dmsp_img = dmsp_img.select(['avg_vis'], ['DMSP'])
        dmsp = ee.ImageCollection([dmsp_img])

        zeros_image = ee.Image(0).reproject(
            ee.Image(dmsp_original.first()).projection()
        )
        zeros_image = zeros_image.set('system:time_start', ee.Date('2014-01-01').millis())
        zeros_image = zeros_image.set('system:time_end', ee.Date('3000-01-01').millis())
        zeros_image_coll = ee.ImageCollection([zeros_image])

        requested_daterange = ee.DateRange(self.start_date, self.end_date)
        viirs_start = ee.Date('2014-1-1')
        use_zeros_image = requested_daterange.contains(viirs_start)

        start_year = ee.Number(ee.Date(self.start_date).get('year'))
        viirs_start_year = ee.Number(viirs_start.get('year'))
        use_zeros_image2 = start_year.gte(viirs_start_year)

        dmsp = ee.Algorithms.If(use_zeros_image, zeros_image_coll, dmsp)
        dmsp = ee.Algorithms.If(use_zeros_image2, zeros_image_coll, dmsp)
        self.dmsp = ee.ImageCollection(dmsp)

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


class ModisDailyLst(GlobalImageDatasource):
    """
    Land Surface Temperature

    Exported Bands:
        LST_DAY
        LST_NIGHT

    All Bands:

    LST_Day_1km	Kelvin	7500	65535	0.02	
    Daytime Land Surface Temperature

    QC_Day					
    Daytime LST Quality Indicators

    Day_view_time	Hours	0	240	0.1	
    Local time of day observation

    Day_view_angle	Degrees	0	130		
    View zenith angle of day observation

    LST_Night_1km	Kelvin	7500	65535	0.02	
    Nighttime Land Surface Temperature

    QC_Night					
    Nighttime LST Quality indicators

    Night_view_time	Hours	0	240	0.1	
    Local time of night observation

    Night_view_angle	Degrees	0	130		
    View zenith angle of night observation

    Emis_31		1	255	0.002	
    Band 31 emissivity

    Emis_32		1	255	0.002	
    Band 32 emissivity

    Clear_day_cov		1	65535	0.0005	
    Day clear-sky coverage

    Clear_night_cov			1	65535	0.0005	
    Night clear-sky coverage
    """
    @staticmethod
    def band_names():
        return ['LST_DAY', 'LST_NIGHT']

    def build_img_coll(self):
        modis = ee.ImageCollection('MODIS/006/MOD11A1')
        self.bands = ['LST_DAY', 'LST_NIGHT']

        self.modis = modis.filterDate(
            self.start_date, self.end_date
        ).map(
            self.qa_filter
        ).select(
            ['LST_Day_1km', 'LST_Night_1km'],
            self.bands
        )

    def get_img_coll(self):
        return self.modis

    @classmethod
    def decode_qa(cls, img):

        def _decode_qa(band_name, prefix):
            qa = img.select(band_name)

            # 0: LST produced, good quality, not necessary to examine more detailed QA
            # 1: LST produced, other quality, recommend examination of more detailed QA
            # 2: LST not produced due to cloud effects
            # 3: LST not produced primarily due to reasons other than cloud
            mandatory_flags = qa.bitwiseAnd(3).select(
                [0], ['{}_MANDATORY_QA'.format(prefix)]
            )
            # 0: Good data quality
            # 1: Other quality data
            # 2: TBD
            # 3: TBD
            data_quality_flags = qa.leftShift(2).bitwiseAnd(3).select(
                [0], ['{}_QUALITY_QA'.format(prefix)]
            )
            # 0: Average emissivity error <= 0.01
            # 1: Average emissivity error <= 0.02
            # 2: Average emissivity error <= 0.04
            # 3: Average emissivity error > 0.04
            emissivity_error_flags = qa.leftShift(4).bitwiseAnd(3).select(
                [0], ['{}_EMISSIVITY_QA'.format(prefix)]
            )
            # 0: Average LST error <= 1K
            # 1: Average LST error <= 2K
            # 2: Average LST error <= 3K
            # 3: Average LST error > 3K
            lst_error_flags = qa.leftShift(6).bitwiseAnd(3).select(
                [0], ['{}_LST_QA'.format(prefix)]
            )

            return ee.Image([
                mandatory_flags,
                data_quality_flags,
                emissivity_error_flags,
                lst_error_flags,
            ])

        day_qa = _decode_qa('QC_Day', 'DAY')
        night_qa = _decode_qa('QC_Night', 'NIGHT')
        return day_qa, night_qa

    @classmethod
    def qa_filter(cls, img):
        day_qa, night_qa = cls.decode_qa(img)
        day_mask = day_qa.select('DAY_MANDATORY_QA').lt(2)
        night_mask = night_qa.select('NIGHT_MANDATORY_QA').lt(2)

        day_lst = img.select('LST_Day_1km').updateMask(day_mask)
        night_lst = img.select('LST_Night_1km').updateMask(night_mask)

        return ee.Image([day_lst, night_lst])


class ModisDailyLstLongRunningAverage(SingleImageDatasource):
    """
    Bands:
        LST_DAY__long_running_average  # Over the life of the collection
        LST_NIGHT__long_running_average
        LST_DAY__month_{}_long_running_average  # Over the life of the collection only on the given month
        LST_NIGHT__month_{}_long_running_average
        + LST_DAY__{}_long_running_average  # For each season in season_to_int_months
        + LST_NIGHT__{}_long_running_average  # For each season in season_to_int_months
    """
    @staticmethod
    def _build_helper(season_to_int_months=None):
        if season_to_int_months is None:
            season_to_int_months = {}

        current_year = datetime.datetime.now().year
        long_running_years = list(range(2000, current_year + 1))

        modis_lst = ModisDailyLst(
            '{}-01-01'.format(min(long_running_years)),
            '{}-12-31'.format(max(long_running_years))
        ).get_img_coll()

        def _process_img_coll(img_coll, suffix):
            original_bands = ['LST_DAY', 'LST_NIGHT']
            bands = [
                '{}__{}'.format(original, suffix)
                for original in original_bands
            ]
            img = img_coll.mean().select(original_bands, bands)
            return img, bands

        long_running_img, bands = compute_monthly_seasonal_annual_and_long_running(
            base_img_coll=modis_lst,
            long_running_years=long_running_years,
            season_to_int_months=season_to_int_months,
            process_img_coll=_process_img_coll,
        )

        return long_running_img, bands

    @staticmethod
    def band_names(season_to_int_months=None):
        _, bands = ModisDailyLstLongRunningAverage._build_helper(
            season_to_int_months=season_to_int_months
        )
        return bands

    def build_img_coll(self, season_to_int_months=None):
        long_running_img, _ = ModisDailyLstLongRunningAverage._build_helper(
            season_to_int_months=season_to_int_months
        )
        self.coll = ee.ImageCollection([long_running_img])

    def get_img_coll(self):
        return self.coll


class ModisNdvi(GlobalImageDatasource):
    """
    Bands:
        NDVI
        EVI
    """
    @staticmethod
    def band_names():
        return ['NDVI', 'EVI']

    def build_img_coll(self):
        ndvi = ee.ImageCollection('MODIS/006/MOD13A2')
        # original_bands = [
        #     'NDVI',
        #     'EVI',
        #     'DetailedQA',
        #     'sur_refl_b01',  # Red
        #     'sur_refl_b02',  # NIR
        #     'sur_refl_b03',  # Blue
        #     'sur_refl_b07',  # mid-IR
        #     'ViewZenith',
        #     'SolarZenith',
        #     'RelativeAzimuth',
        #     'DayOfYear',
        #     'SummaryQA',
        # ]
        self.ndvi = ndvi.filterDate(
            self.start_date, self.end_date
        ).map(
            self.qa_filter
        ).select(
            ['NDVI', 'EVI'],
            ['NDVI', 'EVI']
        )

    def get_img_coll(self):
        return self.ndvi

    @classmethod
    def decode_qa(cls, img):
        qa = img.select('DetailedQA')

        # 0: VI produced with good quality
        # 1: VI produced, but check other QA
        # 2: Pixel produced, but most probably cloudy
        # 3: Pixel not produced due to other reasons than clouds
        vi_quality = qa.bitwiseAnd(1).add(qa.bitwiseAnd(2))
        vi_quality = vi_quality.int().rename(['vi_quality'])

        # 0: Highest quality
        # 1: Lower quality
        # 2: Decreasing quality
        # 4: Decreasing quality
        # 8: Decreasing quality
        # 9: Decreasing quality
        # 10: Decreasing quality
        # 12: Lowest quality
        # 13: Quality so low that it is not useful
        # 14: L1B data faulty
        # 15: Not useful for any other reason/not processed
        vi_usefulness = qa.bitwiseAnd(60).rightShift(ee.Image(2))
        vi_usefulness = vi_usefulness.int().rename(['vi_usefulness'])

        # 0: Climatology
        # 1: Low
        # 2: Intermediate
        # 3: High
        aerosol_quantity = qa.bitwiseAnd(192).rightShift(ee.Image(6))
        aerosol_quantity = aerosol_quantity.int().rename(['aerosol_quantity'])

        adjacent_cloud_detected = qa.bitwiseAnd(256).neq(0)
        adjacent_cloud_detected = adjacent_cloud_detected.updateMask(adjacent_cloud_detected).rename(['adjacent_cloud_detected'])

        atmo_brdf_corr = qa.bitwiseAnd(512).neq(0)
        atmo_brdf_corr = atmo_brdf_corr.updateMask(atmo_brdf_corr).rename(['atmo_brdf_corr'])

        mixed_clouds = qa.bitwiseAnd(1024).neq(0)
        mixed_clouds = mixed_clouds.updateMask(mixed_clouds).rename(['mixed_clouds'])

        # 0: Shallow ocean
        # 1: Land (nothing else but land)
        # 2: Ocean coastlines and lake shorelines
        # 3: Shallow inland water
        # 4: Ephemeral water
        # 5: Deep inland water
        # 6: Moderate or continental ocean
        # 7: Deep ocean
        land_water_mask = qa.bitwiseAnd(14336).rightShift(ee.Image(11))
        land_water_mask = land_water_mask.int().rename(['land_water_mask'])

        poss_snow_ice = qa.bitwiseAnd(16384).neq(0)
        poss_snow_ice = poss_snow_ice.updateMask(poss_snow_ice).rename(['poss_snow_ice'])

        poss_shadow = qa.bitwiseAnd(32768).neq(0)
        poss_shadow = poss_shadow.updateMask(poss_shadow).rename(['poss_shadow'])

        # img.select('SummaryQA')
        # 0		Good Data: use with confidence
        # 1		Marginal Data: useful, but look at other QA information
        # 2		Snow/Ice: target covered with snow/ice
        # 3		Cloudy: target not visible, covered with cloud

        return ee.Image([
            vi_quality, vi_usefulness, aerosol_quantity,
            adjacent_cloud_detected, atmo_brdf_corr,
            mixed_clouds, land_water_mask, poss_snow_ice,
            poss_shadow, img.select('SummaryQA').rename(['summary_qa'])
        ])

    @classmethod
    def qa_filter(cls, img):
        qa = cls.decode_qa(img)
        mask = qa.select('summary_qa').neq(3).And(qa.select('vi_usefulness').lte(12))
        return img.updateMask(mask)
