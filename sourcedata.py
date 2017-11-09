import ee
import lndsatimgtools as imgtools


def _mergejoin(joinedelement):
    #  Inner join returns a FeatureCollection with a primary and secondary set of
    #  properties. Properties are collapsed into different bands of an image.
    return ee.Image.cat(joinedelement.get('primary'), joinedelement.get('secondary'))


#  Convenience function for joining two collections based on system:time_start
def joincoll(coll1, coll2):
    eqfilter = ee.Filter.equals(rightField='system:time_start', leftField='system:time_start')
    join = ee.Join.inner()
    joined = ee.ImageCollection(join.apply(coll1, coll2, eqfilter))
    return joined.map(_mergejoin).sort('system:time_start')


class LandsatTOA:

    def __init__(self, filterpoly, start_date, end_date):

        self.filterpoly = filterpoly

        self.s = start_date
        self.e = end_date

        self.l8 = self._init_coll('LANDSAT/LC8_L1T_TOA_FMASK')
        self.l7 = self._init_coll('LANDSAT/LE7_L1T_TOA_FMASK')
        self.l5 = self._init_coll('LANDSAT/LT5_L1T_TOA_FMASK')

        self.l8fm = self.l8.map(self._mask_fm).map(self._rename_l8)
        self.l7fm = self.l7.map(self._mask_fm).map(self._rename_l7)
        self.l5fm = self.l5.map(self._mask_fm).map(self._rename_l5)

        self.merged = ee.ImageCollection(self.l5.merge(self.l7).merge(self.l8)).sort('system:time_start')
        self.mergedfm = ee.ImageCollection(self.l5fm.merge(self.l7fm).merge(self.l8fm)).sort('system:time_start')

    @staticmethod
    def _pansharpen(image):

        # Convert the RGB bands to the HSV color space.
        rgb = image.select('RED', 'BLUE', 'GREEN')
        gray = image.select('PAN')

        # Swap in the panchromatic band and convert back to RGB.
        hsv = rgb.rgbToHsv().select('hue', 'saturation')

        rgb_sharpened = ee.Image.cat(hsv, gray).hsvToRgb()

        return rgb_sharpened

    def _init_coll(self, name):

        return ee.ImageCollection(name).filterBounds(self.filterpoly).filterDate(self.s, self.e)

    @staticmethod
    def _mask_fm(img):

        """
        fmask: cloud
        mask
        0 = clear
        1 = water
        2 = shadow
        3 = snow
        4 = cloud
        """

        fmask = img.select('fmask')
        # goodpx = fmask.lt(2)
        goodpx = fmask.eq(0)

        return img.updateMask(goodpx).select(img.bandNames().removeAll(['fmask', 'BQA']))

    @staticmethod
    def _rescale_all_bands(img, reflbands, thermbands):

        reflimg = img.select(reflbands).multiply(1000).toInt16()
        thermimg = img.select(thermbands).multiply(10).toInt16()

        return reflimg.addBands(thermimg)

    @staticmethod
    def _rename_l8(l8img):

        """
        Bands of Landsat 8 are:
        B1: Coastal aerosol (0.43 - 0.45 um)
        B2: Blue (0.45 - 0.51 um)
        B3: Green (0.53 - 0.59 um)
        B4: Red (0.64 - 0.67 um)
        B5: Near Infrared (0.85 - 0.88 um)
        B6: Short-wave Infrared 1 (1.57 - 1.65 um)
        B7: Short-wave infrared 2 (2.11 - 2.29 um)
        B8: Panchromatic (0.50 - 0.68 um)
        B9: Cirrus (1.36 - 1.38 um)
        B10: Thermal Infrared 1 (10.60 - 11.19 um)
        B11: Thermal Infrared 2 (11.50 - 12.51 um)
        :param l8img:
        :return:
        """

        bands = ['AEROS', 'BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'PAN', 'CIRR', 'THER1', 'THER2']
        reflimg = l8img.select([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], bands)

        return reflimg

    @staticmethod
    def _rename_l7(l7img):

        """
        B1 - blue	0.45 - 0.52
        B2 - green	0.52 - 0.60
        B3 - red	0.63 - 0.69
        B4 - Near Infrared	0.77 - 0.90
        B5 - Short-wave Infrared	1.55 - 1.75
        B6_VCID_1 - Thermal Infrared	10.40 - 12.50
        B6_VCID_2 - Thermal Infrared	10.40 - 12.50
        B7 - Short-wave Infrared	2.09 - 2.35
        B8 - Panchromatic (Landsat 7 only)	0.52 - 0.90
        """

        reflbands = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'PAN']
        reflimg = l7img.select([0, 1, 2, 3, 4, 7, 8], reflbands)

        thermbands = ['THER1', 'THER2']
        thermimg = l7img.select([5, 6], thermbands)

        return reflimg.addBands(thermimg)

    @staticmethod
    def _rename_l5(l5img):

        """
        B1 - blue	0.45 - 0.52
        B2 - green	0.52 - 0.60
        B3 - red	0.63 - 0.69
        B4 - Near Infrared	0.77 - 0.90
        B5 - Short-wave Infrared	1.55 - 1.75
        B6 - Thermal Infrared	10.40 - 12.50
        B7 - Short-wave Infrared	2.09 - 2.35
        """

        reflbands = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2']
        reflimg = l5img.select([0, 1, 2, 3, 4, 6], reflbands)

        thermbands = ['THER1']
        thermimg = l5img.select([5], thermbands)

        return reflimg.addBands(thermimg)


class LandsatLEDAPS:
    """
    Landsat surface reflectance images as computed by the LEDAPS method (http://ledaps.nascom.nasa.gov/).
    Reflectance is a unitless ratio rescaled to 0-10000.
    An additional atmos_opacity band is added with a representation of atmospheric opacity due to moisture
    and other factors.
    A QA band is added with the the following indicator bits,
    0:unused,
    1:valid data (0=yes, 1=no),
    2:ACCA cloud bit (1=cloudy, 0=clear),
    3:unused,
    4:ACCA snow mask,
    5:land mask based on DEM (1=land, 0=water),
    6:DDV (Dense Dark Vegetation)
    """
    def __init__(self, filterpoly, start_date, end_date):
        self.filterpoly = filterpoly
        self.s = start_date
        self.e = end_date
        self.l7 = self.init_coll('LEDAPS/LE7_L1T_SR')
        self.l5 = self.init_coll('LEDAPS/LT5_L1T_SR')
        self.l7qam = self.l7.map(self.mask_qa)
        self.l5qam = self.l5.map(self.mask_qa)
        self.merged = ee.ImageCollection(self.l5.merge(self.l7)).sort('system:time_start')
        self.mergedqam = ee.ImageCollection(self.l5qam.merge(self.l7qam)).sort('system:time_start')

    def init_coll(self, name):
        return ee.ImageCollection(name).filterBounds(self.filterpoly).filterDate(self.s, self.e).map(
            self.rename_l457)

    @staticmethod
    def rename_l457(limg):
        """
        B1 - blue	0.45 - 0.52
        B1 - green	0.52 - 0.60
        B3 - red	0.63 - 0.69
        B4 - Near Infrared	0.77 - 0.90
        B5 - Short-wave Infrared	1.55 - 1.75
        B7 - Short-wave Infrared	2.09 - 2.35
        atmos_opacity - Atmospheric Opacity
        QA - Bit-packed quality masks
        lndcal_QA -
        """
        return limg.rename(['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'AO', 'QA', 'LCQA'])

    @staticmethod
    def mask_qa(ledapsimg):
        """
        The "QA" band has various flags encoded in different bits.  We extract
        some of them as individual mask bands.
        QA Bit 2: Invalid pixel indicator.
        QA Bit 3: Cloud indicator.
        QA Bit 5: Water indicator.  (0 == water).
        QA Bit 6: Pixel used as "dense dark vegetation"
        """
        valid = ledapsimg.select('QA').bitwiseAnd(2).eq(0)
        valid = valid.updateMask(valid)
        nocloud = ledapsimg.select('QA').bitwiseAnd(4).eq(0)
        nocloud = nocloud.updateMask(nocloud)
        # This flag is technically a "not water" flag
        notwater = ledapsimg.select('QA').bitwiseAnd(32).neq(0)
        notwater = notwater.updateMask(notwater)
        dense_dark_vegetation = ledapsimg.select('QA').bitwiseAnd(64).neq(0)
        dense_dark_vegetation = dense_dark_vegetation.updateMask(dense_dark_vegetation)
        totmask = valid.And(nocloud).And(notwater)
        return ledapsimg.updateMask(totmask)


class LandsatSR:
    def __init__(self, filterpoly, start_date, end_date):
        self.filterpoly = filterpoly
        self.s = start_date
        self.e = end_date
        self.l8 = self.init_coll8('LANDSAT/LC8_SR').map(self.fixwrs)
        self.l7 = self.init_coll('LANDSAT/LE7_SR').map(self.fixwrs)
        self.l5 = self.init_coll('LANDSAT/LT5_SR').map(self.fixwrs)
        self.l8cfm = self.l8.map(self.cfmask)
        self.l7cfm = self.l7.map(self.cfmask)
        self.l5cfm = self.l5.map(self.cfmask)
        self.l7qam = self.l7.map(self.mask_qa)
        self.l5qam = self.l5.map(self.mask_qa)
        # todo: update names to actual colors
        self.merged = ee.ImageCollection(self.l5.merge(self.l7).merge(self.l8)).sort('system:time_start')
        self.mergedcfm = ee.ImageCollection(self.l5cfm.merge(self.l7cfm).merge(self.l8cfm)).sort('system:time_start')
        self.mergedqam = ee.ImageCollection(self.l5qam.merge(self.l7qam)).sort('system:time_start')

    def init_coll(self, name):
        return ee.ImageCollection(name).filterBounds(self.filterpoly).filterDate(self.s, self.e).map(self.rename_l457)

    def init_coll8(self, name):
        return ee.ImageCollection(name).filterBounds(self.filterpoly).filterDate(self.s, self.e).map(self.rename_l8)

    @staticmethod
    def fixwrs(srimg):
        return srimg.set({'WRS_PATH': srimg.get('wrs_path'), 'WRS_ROW': srimg.get('wrs_row')})

    @staticmethod
    def rename_l8(l8img):
        """
        Bands of Landsat 8 are:
        B1: Coastal aerosol (0.43 - 0.45 um) (signed int16)
        B2: Blue (0.45 - 0.51 um) (signed int16)
        B3: Green (0.53 - 0.59 um) (signed int16)
        B4: Red (0.64 - 0.67 um) (signed int16)
        B5: Near Infrared (0.85 - 0.88 um) (signed int16)
        B6: Short-wave Infrared 1 (1.57 - 1.65 um) (signed int16)
        B7: Short-wave infrared 2 (2.11 - 2.29 um) (signed int16)
        cfmask - cloud mask (unsigned int8)
        cfmask_conf - cloud mask confidence (unsigned int8)
        :param l8img:
        :return:
        """
        return l8img.rename(['AEROS', 'BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'cfmask', 'cfmask_conf'])

    @staticmethod
    def rename_l457(limg):
        """
        B1 - blue	0.45 - 0.52
        B1 - green	0.52 - 0.60
        B3 - red	0.63 - 0.69
        B4 - Near Infrared	0.77 - 0.90
        B5 - Short-wave Infrared	1.55 - 1.75
        B7 - Short-wave Infrared	2.09 - 2.35
        cfmask - cloud mask
        cfmask_conf - cloud mask confidence
        adjacent_cloud_qa - Binary QA mask (0/255)
        cloud_qa - Binary QA mask (0/255)
        cloud_shadow_qa - Binary QA mask (0/255)
        ddv_qa, fill_qa - Binary QA mask (0/255)
        land_water_qa - Binary QA mask (0/255)
        snow_qa - Binary QA mask (0/255)
        atmos_opacity - Atmospheric Opacity
        """
        return limg.rename(['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'cfmask', 'cfmask_conf',
                            'adjacent_cloud_qa', 'cloud_qa', 'cloud_shadow_qa', 'ddv_qa', 'fill_qa',
                            'land_water_qa', 'snow_qa', 'atmos_opacity'])

    @staticmethod
    def mask_qa(srimg):
        """
        QA bands are only available for Landsat 4-7:
        adjacent_cloud_qa,
        cloud_qa,
        cloud_shadow_qa,
        ddv_qa, (dense dark vegetation)
        fill_qa,
        land_water_qa,
        snow_qa
        255 if the corresponding condition was detected
        0 otherwise
        :param srimg: SR Landsat image
        :return: SR Landsat image masked for clear pixels only
        """
        qabands = ['cloud_qa', 'cloud_shadow_qa', 'adjacent_cloud_qa', 'land_water_qa', 'snow_qa']
        badqamask = srimg.select(qabands).reduce(ee.Reducer.anyNonZero())
        totqamask = ee.Image(1).mask(srimg.select('BLUE').mask()).where(badqamask, 0)
        return srimg.updateMask(totqamask)

    @staticmethod
    def cfmask(srimg):
        """
        cfmask: cloud mask.
        0=clear
        1=water
        2=shadow
        3=snow
        4=cloud
        cfmask_conf: cloud mask confidence
        0=none
        1=cloud confidence >= 12.5
        2=cloud confidence > 12.5% and <= 22.5%
        3=cloud confidence > 22.5
        :param srimg:
        :return:
        """
        cfmask = srimg.select('cfmask')
        clearpx = cfmask.eq(0)
        return srimg.updateMask(clearpx)


class LandsatJoined:
    def __init__(self, filterpoly, start_date, end_date):
        self.toa = LandsatTOA(filterpoly, start_date, end_date)
        self.sr = LandsatSR(filterpoly, start_date, end_date)
        # TODO: does doing three joins separately cost more than a single all-inclusive join?
        self.l8 = joincoll(self.sr.l8sel, self.toa.l8sel).map(self._cast_img2float)
        self.l7 = joincoll(self.sr.l7sel, self.toa.l7sel).map(self._cast_img2float)
        self.l5 = joincoll(self.sr.l5sel, self.toa.l5sel).map(self._cast_img2float)
        self.merged = ee.ImageCollection(self.l5.merge(self.l7).merge(self.l8))
        self.merged = self.merged.sort('system:time_start')

    @staticmethod
    def _cast_img2float(img):
        return img.toFloat()


class LandsatPlus(object):
    """
    **** DEPRECATED: DO NOT USE ****
    """
    def __init__(self, filterpoly, start_date, end_date, fplus, fplusargs):
        self._fplus = fplus
        self.joined = LandsatJoined(filterpoly, start_date, end_date)
        self.filterpoly = filterpoly
        self.l8 = self._fplus(self.joined.l8, **fplusargs)
        self.l7 = self._fplus(self.joined.l7, **fplusargs)
        self.l5 = self._fplus(self.joined.l5, **fplusargs)
        self.merged = ee.ImageCollection(self.l5.merge(self.l7).merge(self.l8))
        self.merged = self.merged.sort('system:time_start')


class LandsatPrecleaned(object):
    """
    **** DEPRECATED: DO NOT USE ****
    """
    def __init__(self, filterpoly, start_date, end_date, maskf, precleanscale=120, filtercloudlevel=0):
        self.maskf = maskf
        self.joined = LandsatJoined(filterpoly, start_date, end_date)
        self.filterpoly = filterpoly
        self.precleanscale = precleanscale
        self.l8 = self.get_plus(self.joined.l8).filter(ee.Filter.gt('CLOUDFREE', filtercloudlevel))
        self.l7 = self.get_plus(self.joined.l7).filter(ee.Filter.gt('CLOUDFREE', filtercloudlevel))
        self.l5 = self.get_plus(self.joined.l5).filter(ee.Filter.gt('CLOUDFREE', filtercloudlevel))
        self.merged = ee.ImageCollection(self.l5.merge(self.l7).merge(self.l8))
        self.merged = self.merged.sort('system:time_start')

    def get_plus(self, coll):
        plus = coll.map(imgtools.add_latlon).map(imgtools.add_gcviband).map(imgtools.add_ndviband)
        return plus.map(self._mask_and_count)

    def _mask_and_count(self, img):
        precount = imgtools.count_region(img, self.filterpoly, self.precleanscale)
        img = img.set({'PRECOUNT': precount})
        img = self.maskf(img)
        postcount = imgtools.count_region(img, self.filterpoly, self.precleanscale)
        img = img.set({
            'POSTCOUNT': postcount,
            'CLOUDFREE': postcount.divide(precount).multiply(100)
        })
        return img


class MODISrefl(object):
    def __init__(self, start_date, end_date):
        self.newnames = ['RED', 'NIR', 'BLUE', 'GREEN', 'SWIR1', 'SWIR2']
        self.orignames = None
        self.collname = None
        self.coll = None
        self.start_date = start_date
        self.end_date = end_date

    def rename(self, img):
        return img.select(self.orignames, self.newnames)


class MODISnbar(MODISrefl):
    """
    MODIS BRDF-adjusted Reflectance 16-day Global 500m
    """
    def __init__(self, start_date, end_date):
        MODISrefl.__init__(self, start_date, end_date)
        self.orignames = ['Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2', 'Nadir_Reflectance_Band3',
                          'Nadir_Reflectance_Band4', 'Nadir_Reflectance_Band6', 'Nadir_Reflectance_Band7']
        self.collname = "MODIS/MCD43A4"
        self.coll = ee.ImageCollection(self.collname).filterDate(start_date, end_date).map(self.rename)


class MODISsr(MODISrefl):
    """
    MODIS Surface Reflectance 8-Day Global 500m
    """
    def __init__(self, start_date, end_date):
        MODISrefl.__init__(self, start_date, end_date)
        self.orignames = ['sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03',
                          'sur_refl_b04', 'sur_refl_b06', 'sur_refl_b07']
        self.collname = "MODIS/MOD09A1"
        self.coll = ee.ImageCollection(self.collname).filterDate(start_date, end_date).map(self.rename)


class Sentinel2TOA(object):

    def __init__(self, filterpoly, start_date, end_date):

        self.filterpoly = filterpoly

        self.s = start_date

        self.e = end_date

    def init_coll(self, name):

        return ee.ImageCollection(name).filterBounds(self.filterpoly).filterDate(self.s, self.e).map(self.rename)

    @staticmethod
    def qa_cloudmask(img):

        # Opaque and cirrus cloud masks cause bits 10 and 11 in QA60 to be set,so values less than 1024 are cloud-free
        mask = ee.Image(0).where(img.select('QA60').gte(1024), 1).Not()

        return img.updateMask(mask)

    @staticmethod
    def _rescale(img, exp, thresholds):

        """
        A helper to apply an expression and linearly rescale the output.
        """

        return img.expression(exp, img=img).subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])

    def add_cloud_score(self, img):

        img = ee.Image(img).divide(1000)

        score = ee.Image(1.0)

        score = score.min(self._rescale(img, 'img.cirrus', [0, 0.1]))

        score = score.min(self._rescale(img, 'img.cb', [0.5, 0.8]))

        score = score.min(self._rescale(img.normalizedDifference(['GREEN', 'SWIR1']), 'img', [0.8, 0.6]))

        # Invert the cloudscore so 1 is least cloudy, and rename the band.
        return img.addBands(ee.Image(1).subtract(score).select([0], ['cloudscore']))

    @staticmethod
    def rename(s2img):

        """
        Band	Use	Wavelength	Resolution
        B1	Aerosols	443nm	60m
        B2	Blue	490nm	10m
        B3	Green	560nm	10m
        B4	Red	665nm	10m
        B5	Red Edge 1	705nm	20m
        B6	Red Edge 2	740nm	20m
        B7	Red Edge 3	783nm	20m
        B8	NIR	842nm	10m
        B8a	Red Edge 4	865nm	20m
        B9	Water vapor	940nm	60m
        B10	Cirrus	1375nm	60m
        B11	SWIR 1	1610nm	20m
        B12	SWIR 2	2190nm	20m
        QA10
        QA20
        QA60
        """

        newnames = ['AEROS', 'BLUE', 'GREEN', 'RED', 'RDED1', 'RDED2', 'RDED3',
                    'NIR', 'RDED4', 'VAPOR', 'CIRRU', 'SWIR1', 'SWIR2', 'QA10',
                    'QA20', 'QA60']

        return s2img.rename(newnames)


class Daymet:
    def __init__(self):
        # NOTE: no filterBounds needed; DAYMET is composed by whole-CONUS images
        self.wholecoll = ee.ImageCollection('NASA/ORNL/DAYMET')

    @staticmethod
    def addsradvp(img):
        """
        Calculate vpd and radiation in units of MJ/m2
        :param img: daymet image
        :return: original daymet image enriched with radn and vpd bands
        """
        sr = img.select('srad')
        dl = img.select('dayl')
        radn = sr.multiply(dl).divide(1000000)
        vpx = img.expression("0.6107 * exp( 17.269*t / (237.3 + t))", {'t': img.select('tmax')})
        vpn = img.expression("0.6107 * exp( 17.269*t / (237.3 + t))", {'t': img.select('tmin')})
        vpd = vpx.subtract(vpn).multiply(0.75)
        img = img.addBands(radn.select([0], ['radn']))
        img = img.addBands(vpd.select([0], ['vpd']))
        return img

    @staticmethod
    def _compute_radn(img):
        sr = img.select('srad')
        dl = img.select('dayl')
        radn = sr.multiply(dl).divide(1000000)
        return radn.select([0], ['radn'])

    @staticmethod
    def _compute_vpd(img):
        vpx = img.expression("0.6107 * exp( 17.269*t / (237.3 + t))", {'t': img.select('tmax')})
        vpn = img.expression("0.6107 * exp( 17.269*t / (237.3 + t))", {'t': img.select('tmin')})
        vpd = vpx.subtract(vpn).multiply(0.75)
        return vpd.select([0], ['vpd'])

    @staticmethod
    def _compute_gdd(img):
        # NOTE: this has a hard-coded base temperature for corn in US.
        gdd_c = img.expression(
            '((30 - (30-Tmx)*(Tmx<30)) + (10 + (Tmn-10)*(Tmn>10)))/2.0 - 10.0',
            {'Tmx': img.select('tmax'), 'Tmn': img.select('tmin')})
        return gdd_c.select([0], ['gddC'])

    def get_mean_radn(self, startdate, enddate):
        c = self.wholecoll.filterDate(ee.Date(startdate), ee.Date(enddate)).map(self._compute_radn)
        return c.mean()

    def get_mean_precip(self, startdate, enddate):
        c = self.wholecoll.filterDate(ee.Date(startdate), ee.Date(enddate))
        return c.select('prcp').mean()

    def get_mean_tmax(self, startdate, enddate):
        c = self.wholecoll.filterDate(ee.Date(startdate), ee.Date(enddate))
        return c.select('tmax').mean()

    def get_mean_vpd(self, startdate, enddate):
        c = self.wholecoll.filterDate(ee.Date(startdate), ee.Date(enddate)).map(self._compute_vpd)
        return c.mean()

    def get_mean_vhinge(self, startdate, enddate):
        vpd = self.get_mean_vpd(startdate, enddate)
        vhinge = vpd.expression("(x-1.6) * (x > 1.6)", {'x': vpd}).select([0], ['vhinge'])
        return vhinge

    def get_mean_phinge(self, startdate, enddate):
        precip = self.get_mean_precip(startdate, enddate)
        phinge = precip.expression("(3-x) * (x < 3)", {'x': precip}).select([0], ['phinge'])
        return phinge

    def get_cumul_gdd(self, startdate, enddate):
        gdd_c = self.wholecoll.filterDate(ee.Date(startdate), ee.Date(enddate)).map(self._compute_gdd)
        gdd_sum_c = gdd_c.sum().select([0], ['gddC'])
        gdd_sum_f = gdd_sum_c.expression('1.8 * x', {'x': gdd_sum_c.select(0)}).select([0], ['gddF'])
        return gdd_sum_c.addBands(gdd_sum_f)

    def get_met_metrics(self, datesdict):
        vpd = self.get_mean_vpd(datesdict['vpd_start'], datesdict['vpd_end'])
        prec = self.get_mean_precip(datesdict['prec_start'], datesdict['prec_end'])
        vhinge = self.get_mean_vhinge(datesdict['vpd_start'], datesdict['vpd_end'])
        phinge = self.get_mean_phinge(datesdict['prec_start'], datesdict['prec_end'])
        radn = self.get_mean_radn(datesdict['radn_start'], datesdict['radn_end'])
        maxt = self.get_mean_tmax(datesdict['tmax_start'], datesdict['tmax_end'])
        gdd_sum = self.get_cumul_gdd(datesdict['gdd_start'], datesdict['gdd_end'])
        return ee.Image.cat(vpd, prec, vhinge, phinge, radn, maxt, gdd_sum)

    def metmetrics_usa(self, yr):
        """
        Calculate monthly and seasonal averages of weather variables (weather metrics).
        These metrics depend on the location and the shape of the yield model trained
        in APSIM.
        :param yr: the year for which weather metrics are computed (int)
        :return: ee.Image() with one band per metric
        """
        # TODO: change this to be mappable (i.e. yr must be a ee.String)
        # yr = ee.String(yr)
        datesdict = dict(vpd_start=ee.Date.fromYMD(yr, 7, 1), vpd_end=ee.Date.fromYMD(yr, 7, 31),
                         prec_start=ee.Date.fromYMD(yr, 6, 1), prec_end=ee.Date.fromYMD(yr, 8, 31),
                         radn_start=ee.Date.fromYMD(yr, 6, 1), radn_end=ee.Date.fromYMD(yr, 8, 31),
                         tmax_start=ee.Date.fromYMD(yr, 8, 1), tmax_end=ee.Date.fromYMD(yr, 8, 31),
                         gdd_start=ee.Date.fromYMD(yr, 4, 1), gdd_end=ee.Date.fromYMD(yr, 10, 15))
        met = self.get_met_metrics(datesdict)
        return met.set({'year': yr})
