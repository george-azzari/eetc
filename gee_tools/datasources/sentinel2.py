import numpy as np
import ee
from gee_tools.datasources.interface import MultiImageDatasource, DatasourceError
from gee_tools.imgtools import appendBand, getGLCMTexture


class Sentinel2TOA(MultiImageDatasource):

    def build_img_coll(self):
        """

        :return:
        """
        self.name = "COPERNICUS/S2"
        self.coll = ee.ImageCollection(self.name).filterBounds(self.filterpoly)
        self.coll = self.coll.filterDate(self.start_date, self.end_date).map(self.rename)

    def get_img_coll(self):

        if self.coll is None:
            raise DatasourceError("Missing collection, make sure name is not None.  Name was {}".format(self.name))
        return self.coll

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

