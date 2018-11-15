import numpy as np
import ee
from gee_tools.datasources.interface import MultiImageDatasource
from gee_tools.imgtools import appendBand, getGLCMTexture


class Sentinel1(MultiImageDatasource):

    def build_img_coll(self):
        """

        :return:
        """
        self.name = 'COPERNICUS/S1_GRD'
        self.coll = ee.ImageCollection(self.name).filterBounds(self.filterpoly)
        self.coll = self.coll.filterDate(self.start_date, self.end_date)

    def get_img_coll(self, correctlia=False, addbands=True, addspeckle=True, addtexture=False, orbit='ascending'):
"""
    Args:
        correctlia (bool):  If True, call self.correctLIA with SRTM terrain features.
            Defaults to False.
        addbands (bool):  If True, add bands for vv - vh and vh / vv with names 'DIFF' and 'RATIO'.  If addspeckle is 
            True, then 'DIFF_RLSPCK' and 'RATIO_RLSPCK' are also added.  Defaults to True.
        addtexture (bool):  TODO.
            Defaults to False.
        orbit (str):  One of 'ascending' or 'descending'.  Filter on the 'orbitProperties_pass' pass property.
            Defaults to 'ascending'.
            
    Returns:
        (ee.ImageCollection):  Sentinel 1 image collection filtered and modified as described by the arguments.
"""
        # Filter to get images from different look angles.
        if orbit == 'ascending':
            orbfilter = ee.Filter.eq('orbitProperties_pass', 'ASCENDING')

        elif orbit == 'descending':
            orbfilter = ee.Filter.eq('orbitProperties_pass', 'DESCENDING')

        else:
          orbfilter = ee.Filter.inList('orbitProperties_pass', ['ASCENDING', 'DESCENDING'])

        iw = self.coll.filter(ee.Filter.And(
          ee.Filter.eq('instrumentMode', 'IW'),
          ee.Filter.eq('resolution', 'H'),
          ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'),
          orbfilter
          )
        )

        if correctlia:
            # Compute terrain features from SRTM and re-project into S1 projection and resolution.
            terrain = ee.Algorithms.Terrain(ee.Image("USGS/SRTMGL1_003"))
            terrain = terrain.select(['elevation', 'slope', 'aspect'], ['ELEV', 'SLO', 'ASP'])
            terrain = terrain.reproject('EPSG:4326', None, 10)
            iw = iw.map(lambda img: self.correctLIA(img, terrain.select('SLO'), terrain.select('ASP')))

        if addspeckle:
            iw = iw.map(lambda img: self.RefinedLeeMulti(img))

        if addbands:
            iw = iw.map(lambda img: self.addBands(img, ''))

            if addspeckle:
                iw = iw.map(lambda img: self.addBands(img, '_RLSPCK'))

        if addtexture:
            iw = iw.map(lambda img: img.addBands(self.getS1Texture(img, 4, '')))

            if addspeckle:
              iw = iw.map(lambda img: img.addBands(self.getS1Texture(img, 4, '_RLSPCK')))

        return iw

    @staticmethod
    def getPolBands(img, suffix=''):
        """
         :param img: Sentinel 1 image
         :return: ee.List with back-scatter bands only

         It only makes sense to convert base bands to natural and back to DB.
         Problem is that sometimes it may be difficult to know a priori what
         polarization bands are available in a given image. This function
         returns a list of the polarization bands actually available for the
         input image.

         This is based on removeAll, which drops elements from the left list
         IF they are in the right list, otherwise they're ignored. So, by passing
         all possible polarization bands to the right, I get all "bad" bands. In
         turn, I can pass that to drop all the bad ones and keep the good ones that
         are in the image.

         """
        # TODO: this function returns actually back-scatter bands; should consider renaming.

        # Get list of band names available in the image
        allbands = img.bandNames()

        # List standard bands (account for suffixes (e.g. in case this is a composite image)
        possiblebands = ee.List(['VV', 'VH', 'HH', 'HV', 'DIFF'])
        possiblebands = possiblebands.map(lambda s: ee.String(s).cat(ee.String(suffix)))

        # List standard speckle-corrected bands (account for suffixes (e.g. in case this is a composite image)
        spcklbands = ee.List(['VV_RLSPCK', 'VH_RLSPCK', 'HH_RLSPCK', 'HV_RLSPCK', 'DIFF_RLSPCK'])
        spcklbands = spcklbands.map(lambda s: ee.String(s).cat(ee.String(suffix)))

        # Remove band names in source image that don't belong to standard set.
        badbands = allbands.removeAll(possiblebands.cat(spcklbands))

        # Return only band names that belong to standard set.
        return allbands.removeAll(badbands)

    @staticmethod
    def toNatural(img):
        """
        Convert dB backscatter image into natural values.
        :param img:
        :return:
        """

        # Note: dummy is used to keep the original footprint (as opposed to start directly with an ee.Image(10)).
        dummy = img.select(0).multiply(0).add(10)

        # For some reason properties are lost without copyProperties.
        return ee.Image(dummy.pow(img.divide(10)).copyProperties(img))

    @staticmethod
    def toDB(img):
        """
        Convert backscatter image to dB.
        :param img:
        :return:
        """

        return ee.Image(img).log10().multiply(10.0)

    def replace_db_bands(self, image, suffix=''):
        # Get the back-scatter bands (accounting for suffix, if any).
        polbands = self.getPolBands(image, suffix)

        # Convert back-scatter bands to natural units.
        return image.addBands(self.toNatural(image.select(polbands)), None, True)

    def getS1Texture(self, image, radius, bandsufx=''):
        """
        NOTE: user should make sure that backscatter bands are in NATURAL units (as opposed to dB)
        :param image:
        :param radius:
        :param extrabands:
        :param bandsufx:
        :return:
        """

        # Convert back-scatter bands to natural units.
        image = self.replace_db_bands(image, bandsufx)

        # Get the back-scatter bands (accounting for suffix, if any).
        polbands = self.getPolBands(image, bandsufx)

        # Define scaling dictionary (GLCM wants int images).
        scaler = ee.Dictionary.fromLists(polbands, ee.List.repeat(1e4, polbands.size()))

        # Get GLCM bands, copy properties from source and return (setting to default kernel and average mode).
        glcmtxt = getGLCMTexture(image, radius, None, True, scaler)

        return ee.Image(glcmtxt.copyProperties(image))

    @staticmethod
    def addLocalViewAngle(image, slope, aspect):
        """
        Compute local view angle.
        Will fill in Google docstring format later.
        :param image:
        :param slope:
        :param aspect:
        :return:
        """

        # We can use the gradient of the "angle" band of the S1 image to derive the S1 azimuth angle.
        inc = image.select('angle')
        azimuth = ee.Terrain.aspect(inc).reduceRegion(ee.Reducer.mean(), inc.get('system:footprint'), 1000)
        azimuth = ee.Number(azimuth.get('aspect'))

        # Here we derive the terrain slope and aspect, and then the projection of the slope.
        slope_projected = slope.multiply(ee.Image.constant(azimuth).subtract(aspect).multiply(np.pi/180).cos())

        # And finally the local incidence angle
        lia = inc.subtract(ee.Image.constant(90).subtract(ee.Image.constant(90).subtract(slope_projected))).abs()

        return image.addBands(lia.select([0], ['LIA']))

    def correctLIA(self, image, slope, aspect):
        """
        Apply Local Incidence Angle corrections.
        Will fill in Google docstring format later.
        :param image:
        :param slope:
        :param aspect:
        :return:
        """

        image = self.addLocalViewAngle(image, slope, aspect)
        # vh = exports.toNatural(image.select(ee.String('VH')))
        # vv = exports.toNatural(image.select(ee.String('VV')))
        lia = image.select(ee.String('LIA'))

        vh = image.select('VH').subtract(lia.multiply(np.pi/180.0).cos().log10().multiply(10.0)).select([0],['VH'])
        vv = image.select('VV').subtract(lia.multiply(np.pi/180.0).cos().log10().multiply(10.0)).select([0],['VV'])

        return image.addBands(vh, None, True).addBands(vv, None, True)

    def addBands(self, img, sufx=''):
        """
        Add extra bands, such as difference and ration of VV and VH.
        Will fill in Google docstring format later.
        :param img:
        :param sufx:
        :return:
        """

        # Get basic backscatter bands (accounting for possible suffixes)
        vh = self.toNatural(img.select(ee.String('VH').cat(ee.String(sufx))))
        vv = self.toNatural(img.select(ee.String('VV').cat(ee.String(sufx))))

        # Compute polarizations difference and ratio, and back to source image.
        dpol = self.toDB(vv.subtract(vh).select([0], [ee.String('DIFF').cat(ee.String(sufx))]))
        fpol = vh.divide(vv).select([0], [ee.String('RATIO').cat(ee.String(sufx))])

        return img.addBands(dpol).addBands(fpol)

    @staticmethod
    def _RefinedLee(natimg):
        """
        The RL speckle filter * for a SINGLE BAND *
        NOTE: natimg must be in natural units, i.e. not in dB!
        :param natimg:
        :return:
        """
        # Set up 3x3 kernels
        weights3 = ee.List.repeat(ee.List.repeat(1,3),3)
        kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False)

        mean3 = natimg.reduceNeighborhood(ee.Reducer.mean(), kernel3)
        variance3 = natimg.reduceNeighborhood(ee.Reducer.variance(), kernel3)

        # Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
        sample_weights = ee.List([
            [0,0,0,0,0,0,0],
            [0,1,0,1,0,1,0],
            [0,0,0,0,0,0,0],
            [0,1,0,1,0,1,0],
            [0,0,0,0,0,0,0],
            [0,1,0,1,0,1,0],
            [0,0,0,0,0,0,0]
        ])

        sample_kernel = ee.Kernel.fixed(7, 7, sample_weights, 3, 3, False)

        # Calculate mean and variance for the sampled windows and store as 9 bands
        sample_mean = mean3.neighborhoodToBands(sample_kernel)
        sample_var = variance3.neighborhoodToBands(sample_kernel)

        # Determine the 4 gradients for the sampled windows
        gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
        gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs())
        gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs())
        gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs())

        # And find the maximum gradient amongst gradient bands
        max_gradient = gradients.reduce(ee.Reducer.max())

        # Create a mask for band pixels that are the maximum gradient
        gradmask = gradients.eq(max_gradient)

        # duplicate gradmask bands: each gradient represents 2 directions
        gradmask = gradmask.addBands(gradmask)

        # Determine the 8 directions
        # TODO: need to re-format style here.
        directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
        directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
        directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
        directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))

        # The next 4 are the not() of the previous 4
        directions = directions.addBands(directions.select(0).Not().multiply(5))
        directions = directions.addBands(directions.select(1).Not().multiply(6))
        directions = directions.addBands(directions.select(2).Not().multiply(7))
        directions = directions.addBands(directions.select(3).Not().multiply(8))

        # Mask all values that are not 1-8
        directions = directions.updateMask(gradmask)

        # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's
        # directional band, and is otherwise masked)
        directions = directions.reduce(ee.Reducer.sum())

        sample_stats = sample_var.divide(sample_mean.multiply(sample_mean))

        # Calculate localNoiseVariance
        sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0])

        # Set up the 7*7 kernels for directional statistics
        rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4))

        diag_weights = ee.List([
            [1,0,0,0,0,0,0],
            [1,1,0,0,0,0,0],
            [1,1,1,0,0,0,0],
            [1,1,1,1,0,0,0],
            [1,1,1,1,1,0,0],
            [1,1,1,1,1,1,0],
            [1,1,1,1,1,1,1]
        ])

        rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, False)
        diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, False)

        # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
        dir_mean = natimg.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1))
        dir_var = natimg.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1))

        dir_mean = dir_mean.addBands(natimg.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)))
        dir_var = dir_var.addBands(natimg.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)))

        # and add the bands for rotated kernels
        # TODO: why is this for-loop here in the original JS implementation? For-loops in mapped functions are not cool.
        for i in [1, 2, 3, 4]:
            dir_mean = dir_mean.addBands(natimg.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
            dir_var = dir_var.addBands(natimg.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
            dir_mean = dir_mean.addBands(natimg.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
            dir_var = dir_var.addBands(natimg.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))

        # "collapse" the stack into a single band image (due to masking, each pixel has just one value in
        # it's directional band, and is otherwise masked)
        dir_mean = dir_mean.reduce(ee.Reducer.sum())
        dir_var = dir_var.reduce(ee.Reducer.sum())

        # A finally generate the filtered value
        varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))
        b = varX.divide(dir_var)
        result = dir_mean.add(b.multiply(natimg.subtract(dir_mean)))

        return result.arrayFlatten([['sum']])

    def _apply_rl(self, s1dbimg, band):
        """
        Convenience function to map RL speckle filter.
        :param s1dbimg: full s1 image in dB units.
        :param band: band to apply speckle filter on
        :return:
        """

        band = ee.String(band)
        s1natimg = self.toNatural(ee.Image(s1dbimg).select(band))
        spck = self._RefinedLee(s1natimg)

        return self.toDB(spck.select([0], [band.cat('_RLSPCK')]))

    def RefinedLeeMulti(self, dbimg):
        """

        :param dbimg:
        :return:
        """

        polbands = self.getPolBands(dbimg, '')

        splist = polbands.map(lambda b: self._apply_rl(dbimg, b))
        spimg = ee.Image(ee.ImageCollection.fromImages(splist).iterate(appendBand))

        return dbimg.addBands(spimg)



