import numpy as np
import ee
from gee_tools.datasources.interface import MultiImageDatasource, DatasourceError
from gee_tools.imgtools import appendBand, getGLCMTexture, addDOY


class Sentinel2TOA(MultiImageDatasource):

    def build_img_coll(self):
        """

        :return:
        """
        self.name = "COPERNICUS/S2"
        self.coll = ee.ImageCollection(self.name).filterBounds(self.filterpoly)
        self.coll = self.coll.filterDate(self.start_date, self.end_date).map(self.rename)

    def get_img_coll(self, addVIs, addRDVIs, addCloudMasks):

        s2 = self.coll

        # .map(tools.addDOY); # NOTE: DOY is necesseary if RF-based cloud mask has to be applied.

        if addCloudMasks:
            s2 = s2.map(self.addAllQAMaps)

        if addVIs:
            s2 = s2.map(self.addSWVIs)

        if addRDVIs:
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
        qa60 = img.select('QA60').updateMask(img.select('CIRRU'))

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

        # Remap the native QA band to same class numbers than the Hollstein and FSE ones.
        img = self.decodeQA60(img)
        img = self.remapQA60(img, 1, 4, 5)

        # Get the Hollstein band.
        hls = self.applyHollsteinTree(img).select([0], ['QA_HOLLST'])

        # Get FSE (version 1, hard coded for now).
        # (our tree needs the DOY band)
        img = addDOY(img)
        fse = self.applyFSETree_V1(img).select([0], ['QA_FSEV1'])

        return img.addBands(hls).addBands(fse)

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
        Name  Min Max   Scale   Resolution	Wavelength	Description

        B1    0   10000	0.0001  60 METERS   443nm       Aerosols

        B2    0   10000	0.0001	10 METERS   490nm       Blue

        B3    0   10000	0.0001	10 METERS   560nm       Green

        B4    0   10000	0.0001	10 METERS   665nm       Red

        B5    0   10000	0.0001	20 METERS   705nm       Red Edge 1

        B6    0   10000	0.0001	20 METERS   740nm       Red Edge 2

        B7    0   10000	0.0001	20 METERS   783nm       Red Edge 3

        B8    0   10000	0.0001	10 METERS   842nm       NIR

        B8a   0   10000	0.0001	20 METERS   865nm       Red Edge 4

        B9    0   10000	0.0001	60 METERS   940nm       Water vapor

        B10   0   10000	0.0001	60 METERS   1375nm      Cirrus

        B11   0   10000	0.0001	20 METERS   1610nm      SWIR 1

        B12   0   10000	0.0001	20 METERS   2190nm      SWIR 2

        QA10                    10 METERS               Always empty

        QA20                    20 METERS               Always empty

        QA60                    60 METERS               Cloud mask
        """

        newnames = ['AEROS', 'BLUE', 'GREEN', 'RED', 'RDED1', 'RDED2', 'RDED3',
                    'NIR', 'RDED4', 'VAPOR', 'CIRRU', 'SWIR1', 'SWIR2', 'QA10',
                    'QA20', 'QA60']

        return s2img.rename(newnames)

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

        return img.addBands(rdndvi1).addBands(rdndvi2).addBands(ndvi)

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

    @staticmethod
    def _diff(img, b1, b2):
        diffb1_b2 = img.select(b1).subtract(img.select(b2))
        return diffb1_b2

    @staticmethod
    def _ratio(img, b1, b2):
        ratiob1_b2 = img.select(b1).divide(img.select(b2))
        return ratiob1_b2

    def applyHollsteinTree(self, img):
        """
        ------------------------------------------------
        Other authors contributed to this method:
            Stefania Di Tommaso (sdom@stanford.edu)
            Calum You (zedseayou@gmail.com)
        Center on Food Security and the Environment
        Department of Earth System Science
        Stanford University

        Derived from:
        Hollstein, A., K. Segl, L. Guanter, M. Brell, and M. Enesco (2016)
        Ready-to-Use Methods for the Detection of Clouds, Cirrus, Snow, Shadow,
        Water and Clear Sky Pixels in Sentinel-2 MSI Images.
        Remote Sensing, 8(8), 666-18, doi:10.3390/rs8080666.
        ------------------------------------------------
        """

        pred = img.select('NIR').multiply(0)

        clear1 = pred.where(img.select('GREEN').lt(0.319*10000),
            (pred.where(img.select('RDED4').lt(0.166*10000),
                (pred.where(self._diff(img,'GREEN','RDED3').lt(0.027*10000),
                    (pred.where(self._diff(img,'VAPOR','SWIR1').lt(-0.097*10000),1))# 0.097 has a minus in front
                )))))

        shadow1 = pred.where(img.select('GREEN').lt(0.319*10000),
            (pred.where(img.select('RDED4').lt(0.166*10000),
                (pred.where(self._diff(img,'GREEN','RDED3').lt(0.027*10000),
                    (pred.where(self._diff(img,'VAPOR','SWIR1').gte(-0.097*10000),2))# 0.097 has a minus in front
                )))))

        water = pred.where(img.select('GREEN').lt(0.319*10000),
            (pred.where(img.select('RDED4').lt(0.166*10000),
                (pred.where(self._diff(img,'GREEN','RDED3').gte(0.027*10000),
                    (pred.where(self._diff(img,'VAPOR','SWIR1').lt(0.021*10000),3))
                )))))

        shadow2 = pred.where(img.select('GREEN').lt(0.319*10000),
            (pred.where(img.select('RDED4').lt(0.166*10000),
                (pred.where(self._diff(img,'GREEN','RDED3').gte(0.027*10000),
                    (pred.where(self._diff(img,'VAPOR','SWIR1').gte(0.021*10000),2))
                )))))

        clear2 = pred.where(img.select('GREEN').lt(0.319*10000),
            (pred.where(img.select('RDED4').gte(0.166*10000),
                (pred.where(self._ratio(img,'BLUE','CIRRU').lt(14.689),
                    (pred.where(self._ratio(img,'BLUE','VAPOR').lt(0.788),1))
                )))))

        cirrus1 = pred.where(img.select('GREEN').lt(0.319*10000),
            (pred.where(img.select('RDED4').gte(0.166*10000),
                (pred.where(self._ratio(img,'BLUE','CIRRU').lt(14.689),
                    (pred.where(self._ratio(img,'BLUE','VAPOR').gte(0.788),4))
                )))))

        clear3 = pred.where(img.select('GREEN').lt(0.319*10000),
            (pred.where(img.select('RDED4').gte(0.166*10000),
                (pred.where(self._ratio(img,'BLUE','CIRRU').gte(14.689),1))
            )))

        cloud1 = pred.where(img.select('GREEN').gte(0.319*10000),
            (pred.where(self._ratio(img,'RDED1','SWIR1').lt(4.330),
                (pred.where(self._diff(img,'SWIR1','CIRRU').lt(0.255*10000),
                    (pred.where(self._diff(img,'RDED2','RDED3').lt(-0.016*10000),5))  # 0.016 has a minus in front
                )))))

        cirrus2 = pred.where(img.select('GREEN').gte(0.319*10000),
            (pred.where(self._ratio(img,'RDED1','SWIR1').lt(4.330),
                (pred.where(self._diff(img,'SWIR1','CIRRU').lt(0.255*10000),
                    (pred.where(self._diff(img,'RDED2','RDED3').gte(-0.016*10000),4))  # 0.016 has a minus in front
                )))))

        clear4 = pred.where(img.select('GREEN').gte(0.319*10000),
            (pred.where(self._ratio(img,'RDED1','SWIR1').lt(4.330),
                (pred.where(self._diff(img,'SWIR1','CIRRU').gte(0.255*10000),
                    (pred.where(img.select('AEROS').lt(0.3*10000),1))
                )))))

        cloud2 = pred.where(img.select('GREEN').gte(0.319*10000),
            (pred.where(self._ratio(img,'RDED1','SWIR1').lt(4.330),
                (pred.where(self._diff(img,'SWIR1','CIRRU').gte(0.255*10000),
                    (pred.where(img.select('AEROS').gte(0.3*10000),5))
                )))))

        clear5 = pred.where(img.select('GREEN').gte(0.319*10000),
            (pred.where(self._ratio(img,'RDED1','SWIR1').gte(4.330),
                (pred.where(img.select('GREEN').lt(0.525*10000),
                    (pred.where(self._ratio(img,'AEROS','RDED1').lt(1.184),1))
                )))))

        shadow3 = pred.where(img.select('GREEN').gte(0.319*10000),
            (pred.where(self._ratio(img,'RDED1','SWIR1').gte(4.330),
                (pred.where(img.select('GREEN').lt(0.525*10000),
                    (pred.where(self._ratio(img,'AEROS','RDED1').gte(1.184),2))
                )))))

        snow = pred.where(img.select('GREEN').gte(0.319*10000),
            (pred.where(self._ratio(img,'RDED1','SWIR1').gte(4.330),
                (pred.where(img.select('GREEN').gte(0.525*10000),6))
            )))

        clear = clear1.Or(clear2).Or(clear3).Or(clear4).Or(clear5)
        shadow = shadow1.Or(shadow2).Or(shadow3)
        # water
        cirrus = cirrus1.Or(cirrus2)
        cloud = cloud1.Or(cloud2)
        # snow

        # TODO: a bit of an ugly line. Parenthesis are critical here, careful when re-formatting.
        final = (((((pred.where(clear.gte(1),1)).where(shadow.gte(1),2)).where(water.gte(1),3)).where(cirrus.gte(1),4)).where(cloud.gte(1),5)).where(snow.gte(1),6)

        return final

    @staticmethod
    def applyFSETree_V1(img):
        """
        Other authors contributed to this method:
            - Stefania Di Tommaso (sdom@stanford.edu)
            - Calum You (zedseayou@gmail.com)
        Center on Food Security and the Environment
        Department of Earth System Science
        Stanford University

        #decion tree built from R text code
        #clear 1
        #shadow 2
        #water 3
        #cirrus 4
        #cloud 5
        #snow 6

        #clear 1
        #unclear 2

        :param img:
        :return:
        """

        treeList = [ # bandnames->capital & classes->numbers
            '1) root 20219 15046 4 (0.2461546071 0.2504080320 0.2475889015 0.2558484594)',
            '     2) AEROS>=2631.5 4561   165 5 (0.0006577505 0.9638237229 0.0030695023 0.0324490243)',
            '       4) AEROS>=3011 4095    35 5 (0.0000000000 0.9914529915 0.0021978022 0.0063492063) *',
            '       5) AEROS< 3011 466   130 5 (0.0064377682 0.7210300429 0.0107296137 0.2618025751)',
            '        10) VAPOR>=1022.5 320    50 5 (0.0000000000 0.8437500000 0.0062500000 0.1500000000) *',
            '        11) VAPOR< 1022.5 146    72 4 (0.0205479452 0.4520547945 0.0205479452 0.5068493151)',
            '          22) DOY< 225 87    24 5 (0.0000000000 0.7241379310 0.0344827586 0.2413793103) *',
            '          23) DOY>=225 59     6 4 (0.0508474576 0.0508474576 0.0000000000 0.8983050847) *',
            '     3) AEROS< 2631.5 15658 10633 4 (0.3176650913 0.0425980330 0.3188146634 0.3209222123)',
            '      6) SWIR1>=1481.5 10594  6038 1 (0.4300547480 0.0608835190 0.0919388333 0.4171228998)',
            '        12) AEROS< 1561.5 5060  1591 1 (0.6855731225 0.0065217391 0.1292490119 0.1786561265)',
            '          24) CIRRU< 153.5 4738  1308 1 (0.7239341494 0.0069649641 0.1266357113 0.1424651752)',
            '            48) VAPOR>=325.5 4091   957 1 (0.7660718651 0.0075776094 0.0652652163 0.1610853092)',
            '              96) AEROS< 1369.5 2236   241 1 (0.8922182469 0.0026833631 0.0402504472 0.0648479428) *',
            '              97) AEROS>=1369.5 1855   716 1 (0.6140161725 0.0134770889 0.0954177898 0.2770889488)',
            '               194) CIRRU< 20.5 1076   288 1 (0.7323420074 0.0046468401 0.0929368030 0.1700743494)',
            '                388) DOY>=240.5 406    33 1 (0.9187192118 0.0000000000 0.0197044335 0.0615763547) *',
            '                 389) DOY< 240.5 670   255 1 (0.6194029851 0.0074626866 0.1373134328 0.2358208955)',
            '                   778) SWIR1>=2618.5 457   108 1 (0.7636761488 0.0087527352 0.0525164114 0.1750547046) *',
            '                   779) SWIR1< 2618.5 213   135 4 (0.3098591549 0.0046948357 0.3192488263 0.3661971831)',
            '                    1558) VAPOR< 482.5 130    64 2 (0.3000000000 0.0000000000 0.5076923077 0.1923076923) *',
            '                    1559) VAPOR>=482.5 83    30 4 (0.3253012048 0.0120481928 0.0240963855 0.6385542169) *',
            '               195) CIRRU>=20.5 779   428 1 (0.4505776637 0.0256739409 0.0988446727 0.4249037227)',
            '                 390) DOY>=176.5 587   287 1 (0.5110732538 0.0306643952 0.1192504259 0.3390119250)',
            '                   780) VAPOR< 438 158    69 1 (0.5632911392 0.0316455696 0.2974683544 0.1075949367) *',
            '                   781) VAPOR>=438 429   218 1 (0.4918414918 0.0303030303 0.0536130536 0.4242424242)',
            '                    1562) AEROS< 1463.5 199    69 1 (0.6532663317 0.0201005025 0.0351758794 0.2914572864) *',
            '                    1563) AEROS>=1463.5 230   106 4 (0.3521739130 0.0391304348 0.0695652174 0.5391304348)',
            '                      3126) SWIR1>=3370.5 38    10 1 (0.7368421053 0.0000000000 0.0000000000 0.2631578947) *',
            '                      3127) SWIR1< 3370.5 192    78 4 (0.2760416667 0.0468750000 0.0833333333 0.5937500000) *',
            '                 391) DOY< 176.5 192    60 4 (0.2656250000 0.0104166667 0.0364583333 0.6875000000) *',
            '            49) VAPOR< 325.5 647   314 2 (0.4574961360 0.0030911901 0.5146831530 0.0247295209)',
            '              98) DOY>=99 539   248 1 (0.5398886827 0.0000000000 0.4322820037 0.0278293135)',
            '               196) SWIR1>=1846 234    49 1 (0.7905982906 0.0000000000 0.1623931624 0.0470085470) *',
            '               197) SWIR1< 1846 305   110 2 (0.3475409836 0.0000000000 0.6393442623 0.0131147541) *',
            '              99) DOY< 99 108     8 2 (0.0462962963 0.0185185185 0.9259259259 0.0092592593) *',
            '          25) CIRRU>=153.5 322    93 4 (0.1211180124 0.0000000000 0.1677018634 0.7111801242)',
            '            50) DOY>=337.5 59    26 2 (0.4406779661 0.0000000000 0.5593220339 0.0000000000)',
            '             100) CIRRU< 178.5 30     5 1 (0.8333333333 0.0000000000 0.1666666667 0.0000000000) *',
            '             101) CIRRU>=178.5 29     1 2 (0.0344827586 0.0000000000 0.9655172414 0.0000000000) *',
            '            51) DOY< 337.5 263    34 4 (0.0494296578 0.0000000000 0.0798479087 0.8707224335) *',
            '        13) AEROS>=1561.5 5534  2019 4 (0.1964221178 0.1105890857 0.0578243585 0.6351644380)',
            '          26) VAPOR< 495.5 1405   787 1 (0.4398576512 0.0284697509 0.1480427046 0.3836298932)',
            '            52) SWIR1>=2187.5 883   359 1 (0.5934314836 0.0260475651 0.0271800680 0.3533408834)',
            '             104) AEROS< 1760.5 484   110 1 (0.7727272727 0.0165289256 0.0371900826 0.1735537190) *',
            '             105) AEROS>=1760.5 399   171 4 (0.3759398496 0.0375939850 0.0150375940 0.5714285714)',
            '               210) VAPOR< 388 152    58 1 (0.6184210526 0.0263157895 0.0328947368 0.3223684211)',
            '                 420) BLUE>=1708 107    26 1 (0.7570093458 0.0280373832 0.0280373832 0.1869158879) *',
            '                 421) BLUE< 1708 45    16 4 (0.2888888889 0.0222222222 0.0444444444 0.6444444444) *',
            '               211) VAPOR>=388 247    68 4 (0.2267206478 0.0445344130 0.0040485830 0.7246963563) *',
            '            53) SWIR1< 2187.5 522   295 4 (0.1800766284 0.0325670498 0.3524904215 0.4348659004)',
            '             106) VAPOR< 383.5 308   165 2 (0.2467532468 0.0064935065 0.4642857143 0.2824675325)',
            '               212) DOY>=137.5 218    86 2 (0.2339449541 0.0045871560 0.6055045872 0.1559633028) *',
            '               213) DOY< 137.5 90    37 4 (0.2777777778 0.0111111111 0.1222222222 0.5888888889)',
            '                 426) DOY< 103.5 31    10 1 (0.6774193548 0.0322580645 0.2903225806 0.0000000000) *',
            '                 427) DOY>=103.5 59     6 4 (0.0677966102 0.0000000000 0.0338983051 0.8983050847) *',
            '             107) VAPOR>=383.5 214    74 4 (0.0841121495 0.0700934579 0.1915887850 0.6542056075) *',
            '          27) VAPOR>=495.5 4129  1153 4 (0.1135868249 0.1385323323 0.0271252119 0.7207556309)',
            '            54) AEROS>=2150.5 1171   434 4 (0.0204953032 0.3415883860 0.0085397096 0.6293766012)',
            '             108) BLUE>=2541.5 120    32 5 (0.0250000000 0.7333333333 0.0166666667 0.2250000000) *',
            '             109) BLUE< 2541.5 1051   341 4 (0.0199809705 0.2968601332 0.0076117983 0.6755470980)',
            '               218) CIRRU< 370.5 862   334 4 (0.0243619490 0.3538283063 0.0092807425 0.6125290023)',
            '                 436) VAPOR>=870.5 392   156 5 (0.0000000000 0.6020408163 0.0025510204 0.3954081633)',
            '                   872) DOY>=58 321    87 5 (0.0000000000 0.7289719626 0.0031152648 0.2679127726) *',
            '                   873) DOY< 58 71     2 4 (0.0000000000 0.0281690141 0.0000000000 0.9718309859) *',
            '                 437) VAPOR< 870.5 470    97 4 (0.0446808511 0.1468085106 0.0148936170 0.7936170213)',
            '                   874) DOY< 213 205    84 4 (0.0829268293 0.2975609756 0.0292682927 0.5902439024)',
            '                    1748) DOY>=208 38     0 5 (0.0000000000 1.0000000000 0.0000000000 0.0000000000) *',
            '                    1749) DOY< 208 167    46 4 (0.1017964072 0.1377245509 0.0359281437 0.7245508982) *',
            '                   875) DOY>=213 265    13 4 (0.0150943396 0.0301886792 0.0037735849 0.9509433962) *',
            '               219) CIRRU>=370.5 189     7 4 (0.0000000000 0.0370370370 0.0000000000 0.9629629630) *',
            '            55) AEROS< 2150.5 2958   719 4 (0.1504394861 0.0581473969 0.0344827586 0.7569303584)',
            '             110) SWIR1>=3940.5 634   315 4 (0.4589905363 0.0141955836 0.0236593060 0.5031545741)',
            '               220) DOY>=183.5 245    95 1 (0.6122448980 0.0204081633 0.0530612245 0.3142857143)',
            '                 440) SWIR1>=4410.5 119    20 1 (0.8319327731 0.0252100840 0.0336134454 0.1092436975) *',
            '                 441) SWIR1< 4410.5 126    62 4 (0.4047619048 0.0158730159 0.0714285714 0.5079365079)',
            '                   882) AEROS< 1747.5 34     6 1 (0.8235294118 0.0000000000 0.0000000000 0.1764705882) *',
            '                   883) AEROS>=1747.5 92    34 4 (0.2500000000 0.0217391304 0.0978260870 0.6304347826) *',
            '               221) DOY< 183.5 389   147 4 (0.3624678663 0.0102827763 0.0051413882 0.6221079692)',
            '                 442) VAPOR< 686 147    66 1 (0.5510204082 0.0068027211 0.0136054422 0.4285714286)',
            '                   884) SWIR1>=4392.5 64    15 1 (0.7656250000 0.0000000000 0.0000000000 0.2343750000) *',
            '                   885) SWIR1< 4392.5 83    35 4 (0.3855421687 0.0120481928 0.0240963855 0.5783132530)',
            '                    1770) DOY>=99 42    13 1 (0.6904761905 0.0000000000 0.0000000000 0.3095238095) *',
            '                    1771) DOY< 99 41     6 4 (0.0731707317 0.0243902439 0.0487804878 0.8536585366) *',
            '                 443) VAPOR>=686 242    63 4 (0.2479338843 0.0123966942 0.0000000000 0.7396694215) *',
            '             111) SWIR1< 3940.5 2324   404 4 (0.0662650602 0.0701376936 0.0374354561 0.8261617900) *',
            '       7) SWIR1< 1481.5 5064  1046 2 (0.0825434439 0.0043443918 0.7934439179 0.1196682464)',
            '        14) BLUE< 1305.5 4364   615 2 (0.0815765353 0.0025206233 0.8590742438 0.0568285976)',
            '          28) SWIR1>=1204.5 940   343 2 (0.2776595745 0.0031914894 0.6351063830 0.0840425532)',
            '            56) BLUE< 875.5 145    37 1 (0.7448275862 0.0000000000 0.1517241379 0.1034482759) *',
            '            57) BLUE>=875.5 795   220 2 (0.1924528302 0.0037735849 0.7232704403 0.0805031447) *',
            '          29) SWIR1< 1204.5 3424   272 2 (0.0277453271 0.0023364486 0.9205607477 0.0493574766)',
            '            58) CIRRU< 196.5 3352   231 2 (0.0283412888 0.0023866348 0.9310859189 0.0381861575) *',
            '            59) CIRRU>=196.5 72    31 4 (0.0000000000 0.0000000000 0.4305555556 0.5694444444)',
            '             118) SWIR1>=713 32     1 2 (0.0000000000 0.0000000000 0.9687500000 0.0312500000) *',
            '             119) SWIR1< 713 40     0 4 (0.0000000000 0.0000000000 0.0000000000 1.0000000000) *',
            '        15) BLUE>=1305.5 700   342 4 (0.0885714286 0.0157142857 0.3842857143 0.5114285714)',
            '          30) VAPOR< 367 330   133 2 (0.1757575758 0.0090909091 0.5969696970 0.2181818182)',
            '            60) SWIR1< 794.5 70    32 1 (0.5428571429 0.0285714286 0.0428571429 0.3857142857)',
            '             120) CIRRU< 9.5 36     3 1 (0.9166666667 0.0000000000 0.0277777778 0.0555555556) *',
            '             121) CIRRU>=9.5 34     9 4 (0.1470588235 0.0588235294 0.0588235294 0.7352941176) *',
            '            61) SWIR1>=794.5 260    66 2 (0.0769230769 0.0038461538 0.7461538462 0.1730769231) *',
            '          31) VAPOR>=367 370    84 4 (0.0108108108 0.0216216216 0.1945945946 0.7729729730) *'
        ]

        treeString = '\n'.join(treeList)
        classifier = ee.Classifier.decisionTree(treeString)

        classified = img.classify(classifier)

        return classified

    @staticmethod
    def applyFSETree_V2(img):
        """
        Other authors contributed to this method:
            - Stefania Di Tommaso (sdom@stanford.edu)
            - Calum You (zedseayou@gmail.com)
        Center on Food Security and the Environment
        Department of Earth System Science
        Stanford University

        Difference with V1: added more Clear polygon in dark veg area for training (avoid shadows problem).

        #decion tree built from R text code
        #clear 1
        #shadow 2
        #water 3
        #cirrus 4
        #cloud 5
        #snow 6

        #clear 1
        #unclear 2

        :param img:
        :return:
        """

        treeList = [  # bandnames->capital & classes->numbers
            '  1) root 37397 26596 1 (0.2888199588 0.2364628179 0.2338155467 0.2409016766)',
            '    2) AEROS>=2588.5 8420   515 5 (0.0003562945 0.9388361045 0.0024940618 0.0583135392)',
            '      4) AEROS>=2999.5 7448   186 5 (0.0000000000 0.9750268528 0.0014769066 0.0234962406) *',
            '      5) AEROS< 2999.5 972   329 5 (0.0030864198 0.6615226337 0.0102880658 0.3251028807)',
            '       10) VAPOR>=1000.5 660   138 5 (0.0000000000 0.7909090909 0.0045454545 0.2045454545)',
            '         20) DOY>=38.5 602    92 5 (0.0000000000 0.8471760797 0.0049833887 0.1478405316) *',
            '         21) DOY< 38.5 58    12 4 (0.0000000000 0.2068965517 0.0000000000 0.7931034483) *',
            '       11) VAPOR< 1000.5 312   131 4 (0.0096153846 0.3878205128 0.0224358974 0.5801282051)',
            '         22) DOY< 223 151    47 5 (0.0000000000 0.6887417219 0.0331125828 0.2781456954) *',
            '         23) DOY>=223 161    22 4 (0.0186335404 0.1055900621 0.0124223602 0.8633540373) *',
            '    3) AEROS< 2588.5 28977 18179 1 (0.3726403699 0.0323705007 0.3010318528 0.2939572765)',
            '      6) VAPOR>=332.5 20414 11175 1 (0.4525815617 0.0454589987 0.1140393847 0.3879200549)',
            '       12) AEROS< 1525.5 10084  3019 1 (0.7006148354 0.0046608489 0.1666005553 0.1281237604)',
            '         24) SWIR1>=1484.5 7350  1513 1 (0.7941496599 0.0048979592 0.0571428571 0.1438095238)',
            '           48) CIRRU< 35.5 5612   770 1 (0.8627940128 0.0021382751 0.0555951532 0.0794725588) *',
            '           49) CIRRU>=35.5 1738   743 1 (0.5724971231 0.0138089758 0.0621403913 0.3515535098)',
            '             98) AEROS< 1412.5 1065   274 1 (0.7427230047 0.0112676056 0.0535211268 0.1924882629) *',
            '             99) AEROS>=1412.5 673   267 4 (0.3031203566 0.0178306092 0.0757800892 0.6032689450) *',
            '         25) SWIR1< 1484.5 2734  1474 2 (0.4491587418 0.0040234089 0.4608632041 0.0859546452)',
            '           50) DOY< 28 871     0 1 (1.0000000000 0.0000000000 0.0000000000 0.0000000000) *',
            '           51) DOY>=28 1863   603 2 (0.1916264090 0.0059044552 0.6763285024 0.1261406334)',
            '            102) SWIR1>=1233.5 634   347 1 (0.4526813880 0.0126182965 0.3675078864 0.1671924290)',
            '              204) BLUE< 923.5 244    59 1 (0.7581967213 0.0000000000 0.1311475410 0.1106557377) *',
            '              205) BLUE>=923.5 390   189 2 (0.2615384615 0.0205128205 0.5153846154 0.2025641026) *',
            '            103) SWIR1< 1233.5 1229   202 2 (0.0569568755 0.0024410090 0.8356387307 0.1049633849) *',
            '       13) AEROS>=1525.5 10330  3703 4 (0.2104549855 0.0852855760 0.0627299129 0.6415295257)',
            '         26) CIRRU< 21.5 4222  2062 4 (0.4071530081 0.0357650403 0.0454760777 0.5116058740)',
            '           52) AEROS< 1783.5 2138   746 1 (0.6510757717 0.0102899906 0.0579981291 0.2806361085)',
            '            104) SWIR1>=2505 1598   418 1 (0.7384230288 0.0037546934 0.0256570713 0.2321652065) *',
            '            105) SWIR1< 2505 540   311 4 (0.3925925926 0.0296296296 0.1537037037 0.4240740741)',
            '              210) CIRRU< 12.5 199    79 1 (0.6030150754 0.0100502513 0.1356783920 0.2512562814) *',
            '              211) CIRRU>=12.5 341   162 4 (0.2697947214 0.0410557185 0.1642228739 0.5249266862) *',
            '           53) AEROS>=1783.5 2084   524 4 (0.1569097889 0.0619001919 0.0326295585 0.7485604607)',
            '            106) VAPOR< 423.5 344   182 4 (0.4215116279 0.0290697674 0.0784883721 0.4709302326)',
            '              212) SWIR1>=3314.5 89    23 1 (0.7415730337 0.0224719101 0.0000000000 0.2359550562) *',
            '              213) SWIR1< 3314.5 255   114 4 (0.3098039216 0.0313725490 0.1058823529 0.5529411765) *',
            '            107) VAPOR>=423.5 1740   342 4 (0.1045977011 0.0683908046 0.0235632184 0.8034482759)',
            '              214) SWIR1>=4966 60    11 1 (0.8166666667 0.0666666667 0.0166666667 0.1000000000) *',
            '              215) SWIR1< 4966 1680   288 4 (0.0791666667 0.0684523810 0.0238095238 0.8285714286) *',
            '         27) CIRRU>=21.5 6108  1641 4 (0.0744924689 0.1195153897 0.0746561886 0.7313359528)' ,
            '           54) BLUE>=2179.5 712   310 4 (0.0154494382 0.4115168539 0.0084269663 0.5646067416)',
            '            108) DOY>=58 560   271 5 (0.0196428571 0.5160714286 0.0107142857 0.4535714286)',
            '              216) BLUE>=2477 134    29 5 (0.0223880597 0.7835820896 0.0149253731 0.1791044776) *',
            '              217) BLUE< 2477 426   196 4 (0.0187793427 0.4319248826 0.0093896714 0.5399061033) *',
            '            109) DOY< 58 152     4 4 (0.0000000000 0.0263157895 0.0000000000 0.9736842105) *',
            '           55) BLUE< 2179.5 5396  1331 4 (0.0822831727 0.0809859155 0.0833951075 0.7533358043) *',
            '      7) VAPOR< 332.5 8563  2168 2 (0.1820623613 0.0011678150 0.7468177041 0.0699521196)',
            '       14) SWIR1>=1688.5 1270   563 1 (0.5566929134 0.0031496063 0.3401574803 0.1000000000)',
            '         28) DOY>=99 1045   344 1 (0.6708133971 0.0019138756 0.2086124402 0.1186602871) *',
            '         29) DOY< 99 225    11 2 (0.0266666667 0.0088888889 0.9511111111 0.0133333333) *',
            '       15) SWIR1< 1688.5 7293  1330 2 (0.1168243521 0.0008227067 0.8176333470 0.0647195941)',
            '         30) DOY< 28 342     0 1 (1.0000000000 0.0000000000 0.0000000000 0.0000000000) *',
            '         31) DOY>=28 6951   988 2 (0.0733707380 0.0008631852 0.8578621781 0.0679038987)',
            '           62) BLUE< 1264.5 5914   570 2 (0.0629015894 0.0000000000 0.9036185323 0.0334798783)',
            '            124) CIRRU< 102.5 5817   495 2 (0.0629190304 0.0000000000 0.9149045900 0.0221763796) *',
            '            125) CIRRU>=102.5 97    28 4 (0.0618556701 0.0000000000 0.2268041237 0.7113402062) *',
            '           63) BLUE>=1264.5 1037   418 2 (0.1330761813 0.0057859209 0.5969141755 0.2642237223)',
            '            126) SWIR1< 765 236   105 4 (0.2923728814 0.0084745763 0.1440677966 0.5550847458)',
            '              252) VAPOR< 180.5 65    17 1 (0.7384615385 0.0000000000 0.0769230769 0.1846153846) *',
            '              253) VAPOR>=180.5 171    52 4 (0.1228070175 0.0116959064 0.1695906433 0.6959064327) *',
            '            127) SWIR1>=765 801   216 2 (0.0861423221 0.0049937578 0.7303370787 0.1785268414) *'
        ]

        treeString = '\n'.join(treeList)
        classifier = ee.Classifier.decisionTree(treeString)

        classified = img.classify(classifier)

        return classified
