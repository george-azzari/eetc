"""
Author: George Azzari (gazzari@stanford.edu)
Center on Food Security and the Environment
Department of Earth System Science
Stanford University
"""

import ee


class Daymet:
    def __init__(self):
        # NOTE: no filterBounds needed; DAYMET is composed by whole-CONUS images
        self.wholecoll = ee.ImageCollection('NASA/ORNL/DAYMET_V3')

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
