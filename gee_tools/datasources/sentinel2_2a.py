import numpy as np
import ee
from gee_tools.datasources.interface import MultiImageDatasource, DatasourceError
from gee_tools.imgtools import appendBand, getGLCMTexture, addDOY


class Sentinel2SR(MultiImageDatasource):

    def build_img_coll(self, addVIs=True, addCloudMasks=True):
        """
        Args:

        addVIs (bool): If True do two things.
            Use self.refl_scale to convert bands to reflectance units.
            Add the following bands.  'NBR1', 'NBR2', 'STI', 'NDTI', 'CRC', 'REIP', 'ChloropIndx', 'MTCI', 'WDRVI', 'NDVI'
            Defaults to True.
        addCloudMasks (bool):  If True, use addAllQAMaps to add cloud masks.
            Defaults to True.
        """

        self.name = "COPERNICUS/S2_SR"
        self.coll = ee.ImageCollection(self.name).filterBounds(self.filterpoly)
        self.coll = self.coll.filterDate(self.start_date, self.end_date).map(self.rename)

        self.opticalbands = ee.List([
            'AEROS',
            'BLUE',
            'GREEN',
            'RED',
            'RDED1',
            'RDED2',
            'RDED3',
            'NIR',
            'RDED4',
            'VAPOR',
            'SWIR1',
            'SWIR2',
            'AOT',
            'WVP',
            'TCI_R',
            'TCI_G',
            'TCI_B'
        ])

        self._addVIs = addVIs
        self._addCloudMasks = addCloudMasks

    def get_img_coll(self, addVIs=None, addCloudMasks=None):
        """
        Args:
            addVIs (Optional[bool]): If True do two things.
                Use self.refl_scale to convert bands to reflectance units.
                Add the following bands.  'NBR1', 'NBR2', 'STI', 'NDTI', 'CRC', 'REIP', 'ChloropIndx', 'MTCI', 'WDRVI', 'NDVI'
                Defaults to constructor args which defaults to True.
            addCloudMasks (Optional[bool]):  If True, use addAllQAMaps to add cloud masks.
                Defaults to constructor args which defuaults to True.
                
        Returns:
            (ee.ImageCollection):  The sentinel 2 image collection modified by arguments.
        """
        # TODO: I may need to reconsider this method. Two main issues:
        # TODO: VIs should always be computed in reflectance units (i.e. after scaling).
        addVIs = self._addVIs if addVIs is None else addVIs
        addCloudMasks = self._addCloudMasks if addCloudMasks is None else addCloudMasks

        s2 = self.coll

        if addCloudMasks:
            s2 = s2.map(self.addAllQAMaps)

        if addVIs:
            # VIs have to be computed in reflectance units.
            s2 = s2.map(self.refl_scale)
            s2 = s2.map(self.addSWVIs)
            s2 = s2.map(self.addRededgeExtras)

        return s2

    @staticmethod
    def decodeQA60(img):
        """
        Bitmask for QA60.

        Bit 10: Opaque clouds
        0: No opaque clouds
        1: Opaque clouds present

        Bit 11: Cirrus clouds
        0: No cirrus clouds
        1: Cirrus clouds present
        */
        """

        # NOTE: updateMask here is to make sure QA60 has same footprint of image (sometimes it does not).
        qa60 = img.select('QA60').updateMask(img.select('AEROS'))

        # Bit 10: 2^10 = 1024
        cloudBitMask = qa60.bitwiseAnd(ee.Number(2).pow(10).int())
        cloud = cloudBitMask.neq(0).rename(['PXQA60_CLOUD']).toInt()
        cloud = cloud.updateMask(cloud)

        # Bit 11: 2^11 = 2048
        cirrusBitMask = qa60.bitwiseAnd(ee.Number(2).pow(11).int())
        cirrus = cirrusBitMask.neq(0).rename(['PXQA60_CIRRUS']).toInt()
        cirrus = cirrus.updateMask(cirrus)

        clear = cloudBitMask.eq(0).And(cirrusBitMask.eq(0)).rename(['PXQA60_CLEAR']).toInt()
        clear = clear.updateMask(clear)

        return img.addBands([cloud, cirrus, clear])

    def remapQA60(self, img, clearval, cirrusval, cloudval):
        """
        Remap classes in QA60 bands to new specified values. This is useful for comparisons with other QA masks.
        :param img:
        :param clearval:
        :param cirrusval:
        :param cloudval:
        :return:
        """

        decoded = self.decodeQA60(img)
        clear = decoded.select('PXQA60_CLEAR')
        cloud = decoded.select('PXQA60_CLOUD')
        cirrus = decoded.select('PXQA60_CIRRUS')

        # Remapping to common class values. Also, I am going to merge these into a collection,
        # so they need to have same band names.
        cirrus = cirrus.where(cirrus.eq(1), cirrusval).select([0],['QA60_DECODED'])
        cloud = cloud.where(cloud.eq(1), cloudval).select([0],['QA60_DECODED'])
        clear = clear.where(clear.eq(1), clearval).select([0],['QA60_DECODED'])

        qamerged = ee.ImageCollection.fromImages([clear, cirrus, cloud]).sum()

        bnames = img.bandNames().removeAll(['PXQA60_CLEAR', 'PXQA60_CLOUD', 'PXQA60_CIRRUS'])#.add(['QA60_DECODED'])

        classvalues = {'clear': clearval, 'cloud': cloudval, 'cirrus': cirrusval}

        return img.select(bnames).addBands(qamerged).set('classvalues', classvalues)

    def addAllQAMaps(self, img):

        # Remap the native QA band to same class numbers.
        img = self.decodeQA60(img)
        img = self.remapQA60(img, 1, 4, 5)

        return img

    @staticmethod
    def _rescale(img, exp, thresholds):

        """
        A helper to apply an expression and linearly rescale the output.
        """

        return img.expression(exp, img=img).subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])

    def add_cloud_score(self, img):
        """
        Older function from EE team that computes a simple cloud score.
        :param img:
        :return:
        """

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
        Name        Min Max     Scale   Resolution	Wavelength	Description

        B1          0   10000	0.0001  60 METERS   443nm       Aerosols

        B2          0   10000	0.0001	10 METERS   490nm       Blue

        B3          0   10000	0.0001	10 METERS   560nm       Green

        B4          0   10000	0.0001	10 METERS   665nm       Red

        B5          0   10000	0.0001	20 METERS   705nm       Red Edge 1

        B6          0   10000	0.0001	20 METERS   740nm       Red Edge 2

        B7          0   10000	0.0001	20 METERS   783nm       Red Edge 3

        B8          0   10000	0.0001	10 METERS   842nm       NIR

        B8a         0   10000	0.0001	20 METERS   865nm       Red Edge 4

        B9          0   10000	0.0001	60 METERS   940nm       Water vapor

        B11         0   10000	0.0001	20 METERS   1610nm      SWIR 1

        B12         0   10000	0.0001	20 METERS   2190nm      SWIR 2

        AOT         0   10000   0.0001  10 METERS               Aerosol Optical Thickness

        WVP         0   10000   0.0001  10 METERS               Water Vapor Pressure

        SCL         1   11      0.0001  20 METERS               Science Classification Map

        TCI_R       0   10000   0.0001  10 METERS               True Color Image, Red Channel

        TCI_G       0   10000   0.0001  10 METERS               True Color Image, Green Channel

        TCI_B       0   10000   0.0001  10 METERS               True Color Image, Blue Channel

        MSK_CLDPRB  0   100     0.0001  20 METERS               Cloud Probability Map

        MSK_SNWPRB  0   100     0.0001  10 METERS               Snow Probability Map

        QA10                            10 METERS               Always empty

        QA20                            20 METERS               Always empty

        QA60                            60 METERS               Cloud mask
        """

        newnames = ['AEROS', 'BLUE', 'GREEN', 'RED', 'RDED1', 'RDED2', 'RDED3',
                    'NIR', 'RDED4', 'VAPOR', 'SWIR1', 'SWIR2', 'AOT', 'WVP', 'SCL',
                    'TCI_R', 'TCI_G', 'TCI_B', 'MSK_CLDPRB', 'MSK_SNWPRB', 
                    'QA10', 'QA20', 'QA60']

        return s2img.rename(newnames)

    @staticmethod
    def refl_scale(img):
        """
        Scale bands back to original reflectance units.
        TODO: update this method to use generalized scaler in imgtools.
        :param img:
        :return:
        """

        optical10  = ['BLUE',  'GREEN', 'RED', 'NIR', 'AOT', 'WVP', 'TCI_R', 'TCI_G', 'TCI_B']

        optical20  = ['RDED1', 'RDED2', 'RDED3',
                      'RDED4', 'SWIR1', 'SWIR2']

        optical60  = ['AEROS', 'VAPOR']

        scaler = ee.Image.constant(0.0001)

        scaler10 = scaler.updateMask(img.select(optical10[0]))
        scoptical10 = img.select(optical10).multiply(scaler10)

        scaler20 = scaler.updateMask(img.select(optical20[0]))
        scoptical20 = img.select(optical20).multiply(scaler20)

        scaler60 = scaler.updateMask(img.select(optical60[0]))
        scoptical60 = img.select(optical60).multiply(scaler60)

        # Re-add and overwrite.
        img = img.addBands(scoptical10, optical10, True)
        img = img.addBands(scoptical20, optical20, True)
        img = img.addBands(scoptical60, optical60, True)

        return img

    @staticmethod
    def addREIP(image):
        """
        Compute the Red Edge Inflaction Point (REIP).

        :param image:
        :return:
        """

        expr = '700 + 40 * (((RED+RDED3)/2)-RDED1)/(RDED3-RDED1)'

        exprdict = {

            'RED': image.select('RED'),

            'RDED3': image.select('RDED3'),

            'RDED2': image.select('RDED2'),

            'RDED1': image.select('RDED1')
        }

        reip = image.expression(expr, exprdict).select([0], ['REIP'])

        return image.addBands(reip)

    @staticmethod
    def addChloropIndx(img):
        """
        Compute:
            - Rededge Chlorophyll Indexes
            - Green Chlorophyll Index
        :param img:
        :return:
        """

        rdgcvi = img.expression(
            '(nir / reded) - 1',
            {
                'nir': img.select('NIR'),
                'reded': img.select('RDED1')

            }).select([0], ['RDGCVI1'])

        rdgcvi2 = img.expression(
            '(nir / reded) - 1',
            {
                'nir': img.select('NIR'),
                'reded': img.select('RDED2')

            }).select([0], ['RDGCVI2'])

        gcvi = img.expression(
            '(nir / green) - 1',
            {'nir': img.select('NIR'),
             'green': img.select('GREEN')
             }).select([0], ['GCVI'])

        return img.addBands(gcvi).addBands(rdgcvi).addBands(rdgcvi2)

    @staticmethod
    def addMTCI(img):
        """
        Compute the MERIS Terrestrial Chlorophyll Index
        :param img:
        :return:
        """

        mtci = img.expression(
            '(nir - reded) / (reded - red)',
            {
                'nir': img.select('NIR'),
                'reded': img.select('RDED1'),
                'red': img.select('RED')

            }).select([0], ['MTCI'])

        mtci2 = img.expression(
            '(reded2 - reded1) / (reded1 - red)',
            {
                'reded2': img.select('RDED2'),
                'reded1': img.select('RDED1'),
                'red': img.select('RED')

            }).select([0], ['MTCI2'])

        return img.addBands(mtci).addBands(mtci2)

    @staticmethod
    def addWDRVI(img, alpha=0.2):

        wdrvi = img.expression(
            '(a * nir - red) / (red + a * nir) + (1-a)/(1+a)',
            {
             'a': alpha,
             'nir': img.select('NIR'),
             'red': img.select('RED')

            }).select([0], ['WDRVI'])

        gwdrvi = img.expression(
            '(a * nir - green) / (green + a * nir) + (1-a)/(1+a)',
            {
                'a': alpha,
                'nir': img.select('NIR'),
                'green': img.select('GREEN')
            }).select([0], ['GRWDRVI'])

        rewdrvi = img.expression(
            '(a * nir - rded) / (rded + a * nir) + (1-a)/(1+a)',
            {
                'a': alpha,
                'nir': img.select('NIR'),
                'rded': img.select('RDED1')
            }).select([0], ['RDWDRVI'])

        return img.addBands(wdrvi, None, True).addBands(gwdrvi).addBands(rewdrvi)

    @staticmethod
    def addNDVI(img):

        # Classic NDVI
        ndvi = img.expression(
            '(nir - red) / (red + nir)',
            {
                'nir': img.select('NIR'),
                'red': img.select('RED')

            }).select([0], ['NDVI'])

        sndvi = img.expression(
            '(nir - red) / (red + nir + 0.16)',
            {
                'nir': img.select('NIR'),
                'red': img.select('RED')
            }).select([0], ['SNDVI']) 

        # Rededge NDVI 1
        rdndvi1 = img.expression(
            '(nir - reded) / (reded + nir)',
            {
                'nir': img.select('NIR'),
                'reded': img.select('RDED1')

            }).select([0], ['RDNDVI1'])

        # Rededge NDVI 2*/
        rdndvi2 = img.expression(
            '(nir - reded) / (reded + nir)',
            {
                'nir': img.select('NIR'),
                'reded': img.select('RDED2')

            }).select([0], ['RDNDVI2'])

        return img.addBands(rdndvi1).addBands(rdndvi2).addBands(ndvi).addBands(sndvi)

    def addRededgeExtras(self, img, alpha=0.2):

        img = self.addREIP(img)
        img = self.addChloropIndx(img)
        img = self.addMTCI(img)
        img = self.addWDRVI(img, alpha)
        img = self.addNDVI(img)

        return img

    @staticmethod
    def addSWVIs(img):

        nbr1 = img.expression(
            '(nir - swir1) / (nir + swir1)',
            {
                'nir': img.select('NIR'),
                'swir1': img.select('SWIR1')

            }).select([0], ['NBR1'])

        nbr2 = img.expression(
            '(nir - swir2) / (nir + swir2)',
            {
                'nir': img.select('NIR'),
                'swir2': img.select('SWIR2')

            }).select([0], ['NBR2'])

        # Simple tillage index */
        sti = img.expression(
            'swir1/swir2',
            {
                'swir1': img.select('SWIR1'),
                'swir2': img.select('SWIR2')

            }).select([0], ['STI'])

        # NDTI */
        ndti = img.expression(
            '(swir1 - swir2) / (swir1 + swir2)',
            {
                'swir1': img.select('SWIR1'),
                'swir2': img.select('SWIR2')

            }).select([0], ['NDTI'])

        # Modified CRC*/
        crc = img.expression(
            '(swir1 - green) / (swir1 + green)',
            {
                'green': img.select('GREEN'),
                'swir1': img.select('SWIR1')

            }).select([0], ['CRC'])

        return img.addBands(nbr1).addBands(nbr2).addBands(sti).addBands(ndti).addBands(crc)

