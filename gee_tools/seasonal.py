import ee
import arraycomposites
import imgtools


class Collection(object):
    """
    Note: the reason why the main argument is a collection generator rather than a collection itself is that
          I can generate the collection specifically and consistently for the seasonal window and year I am
          interested in. This is to avoid generating large collections for nothing (which becomes a problem
          when the collection is based on a join, such as the case of LandsatPlus). Also, I inlcuded "updating"
          functions so that masking nd filtering based on "seasonal" metadata can be performed internally.
    """
    # Todo: could make class more general and include a year range (for example for DSC-like applications)
    def __init__(self, basecollgen, seasonalyear, startmonth, sdoys, edoys, collkwargs):
        self.seasonalyear = ee.Number(seasonalyear)
        self.startmonth = ee.Number(startmonth)
        self.start_agdoys = ee.List(sdoys)
        self.end_agdoys = ee.List(edoys)
        r = lambda d: ee.Date(self.agricdoy2solardate(d, self.startmonth, self.seasonalyear).get('solardate'))
        self.start_solardates = self.start_agdoys.map(r)
        self.end_solardates = self.end_agdoys.map(r)
        self.basecoll = basecollgen(start_date=ee.Date(self.start_solardates.get(0)),
                                    end_date=ee.Date(self.end_solardates.get(-1)),
                                    **collkwargs)
        self.basecoll = self.basecoll.map(self._set_agricdate)

    def _set_agricdate(self, img):
        solardate = ee.Date(ee.Number(img.get('system:time_start')))
        solardoy = ee.Number(solardate.getRelative('day', 'year')).add(ee.Number(1))
        agricdoy = ee.Number(self.solardate2agricdoy(solardate, self.startmonth).get('agdoy'))
        dummy = img.select(0).multiply(0).long()  # this is a workaround for now (see Noel email)
        timestampimg = dummy.add(img.metadata('system:time_start')).select([0], ['MSTIMESTAMP'])
        solardoyimg = dummy.add(solardoy).select([0], ['SOLARDOY'])
        agricdoyimg = dummy.add(agricdoy).select([0], ['AGRICDOY'])
        img = img.addBands(ee.Image([timestampimg, solardoyimg, agricdoyimg]))
        return ee.Image(img).set({'SOLARDOY': solardoy,
                                  'AGRICDOY': agricdoy,
                                  'SOLARYEAR': ee.Number(solardate.get('year')),
                                  'AGRICYEAR': self.seasonalyear})

    def update_basecoll(self, updater):
        self.basecoll = self.basecoll.map(updater)

    def filter_basecoll(self, thisfilter):
        self.basecoll = self.basecoll.filter(thisfilter)

    @staticmethod
    def agricdoy2solardate(agdoy, startmonth, agyear):
        dy = ee.Number(ee.Algorithms.If(ee.Number(startmonth).gt(ee.Number(6)), -1, 0))
        startdate = ee.Date.fromYMD(ee.Number(agyear).add(dy), startmonth, 1)
        solardate = startdate.advance(ee.Number(agdoy).subtract(1), 'day')
        return ee.Dictionary({'solardate': solardate, 'startdate': startdate})

    @staticmethod
    def solardate2agricdoy(solardate, startmonth):
        solaryear = ee.Number(solardate.get('year'))
        solarm = ee.Number(solardate.get('month'))
        startyear = ee.Number(ee.Algorithms.If(solarm.lt(ee.Number(startmonth)), solaryear.subtract(1), solaryear))
        startdate = ee.Date.fromYMD(startyear, startmonth, 1)
        agdoy = solardate.difference(startdate,'day').add(1).toInt()
        return ee.Dictionary({'agdoy': agdoy, 'startdate':startdate})

    @staticmethod
    def solardate2agricyear(solardate, startmonth):
        # Note: not sure if converting to agyear is that important given my use-cases.
        solaryear = ee.Number(solardate.get('year'))
        solarm = ee.Number(solardate.get('month'))
        yoffset = ee.Number(
            ee.Algorithms.If(ee.Number(startmonth).gt(6),
                             ee.Algorithms.If(solarm.gte(startmonth), 1, 0),
                             ee.Algorithms.If(solarm.lt(startmonth), -1, 0)))
        agyear = solaryear.add(yoffset)
        return agyear

    def filter_by_doy(self, seasonindex):
        fdoy = ee.Filter.And(ee.Filter.gte('AGRICDOY', ee.Number(self.start_agdoys.get(seasonindex))),
                             ee.Filter.lt('AGRICDOY', ee.Number(self.end_agdoys.get(seasonindex))))
        return self.basecoll.filter(fdoy)

    @staticmethod
    def nullify_emptyimages(img):
        final = ee.Algorithms.If(ee.Image(img).bandNames(), ee.Image(img), None)
        return final

    def drop_emptyimages(self, imgcoll):
        return imgcoll.map(self.nullify_emptyimages, True)

    def seasonal_qltymos(self, seasonindex, qltyband):
        """
        """
        seasonalcoll = self.filter_by_doy(seasonindex)
        return ee.Image(seasonalcoll.qualityMosaic(qltyband))

    def seasonal_quantmos(self, seasonindex, qltyband, quantile):
        """
        """
        doycoll = self.filter_by_doy(seasonindex)
        return ee.Image(arraycomposites.nthquantile_mosaic(doycoll, qltyband, quantile))

    def biseasonal_qltymos(self, seasonindex1, seasonindex2, qltyband):
        # TODO: how do I drop empty images in this case (e.g. when only one mosaic is empty)?
        # TODO: ee.Image.cat([None, None]) would return an image with 'constant' and 'constant_1' bands.
        dt1 = self.seasonal_qltymos(seasonindex1, qltyband)
        dt2 = self.seasonal_qltymos(seasonindex2, qltyband)
        props = {'bands_1': ee.List(dt1.bandNames()).length(),
                 'bands_2': ee.List(dt2.bandNames()).length(),
                 'year': self.seasonalyear}
        return ee.Image.cat([imgtools.rename_bands(dt1, '1'), imgtools.rename_bands(dt2, '2')]).set(props)

    def biseasonal_quantmos(self, seasonindex1, seasonindex2, qltyband, quantile):
        # TODO: how do I drop empty images in this case (e.g. when only one mosaic is empty)?
        # TODO: ee.Image.cat([None, None]) would return an image with 'constant' and 'constant_1' bands.
        dt1 = self.seasonal_quantmos(seasonindex1, qltyband, quantile)
        dt2 = self.seasonal_quantmos(seasonindex2, qltyband, quantile)
        props = {'bands_1': ee.List(dt1.bandNames()).length(),
                 'bands_2': ee.List(dt2.bandNames()).length(),
                 'year': self.seasonalyear}
        return ee.Image.cat([imgtools.rename_bands(dt1, '1'), imgtools.rename_bands(dt2, '2')]).set(props)


class WRSGroups(object):
    # TODO: figure out a better way to utilize "updater" (if at all). It's currently called for each WSR subset.
    def __init__(self, basecollgen, wrsboxes, updater, seasonalyear, startmonth, sdoys, edoys, collkwargs):
        self.wrsboxes = wrsboxes
        self.updater = updater
        self.seasargs = dict(basecollgen=basecollgen, seasonalyear=seasonalyear, startmonth=startmonth,
                             sdoys=sdoys, edoys=edoys)
        self.collkwargs = collkwargs

    def _biseasonal_qltymoscoll(self, wrsbox, seasonindex1, seasonindex2, qltyband):
        tmpargs = dict(self.collkwargs, filterpoly=wrsbox.centroid().geometry())
        seascoll = Collection(collkwargs=tmpargs, **self.seasargs)
        seascoll.update_basecoll(self.updater)
        clipper = wrsbox.geometry().buffer(-200)
        mos = seascoll.biseasonal_qltymos(seasonindex1, seasonindex2, qltyband).clip(clipper)
        return mos

    def biseasonal_qltymoscoll(self, seasonindex1, seasonindex2, qltyband):
        _get_mos = lambda box: self._biseasonal_qltymoscoll(ee.Feature(box), seasonindex1, seasonindex2, qltyband)
        moscoll = self.wrsboxes.map(_get_mos)
        filt = ee.Filter.And(ee.Filter.neq('bands_1', 0), ee.Filter.neq('bands_2', 0))
        # dropping composites without full set of bands
        return moscoll.filter(filt)

    def _biseasonal_quantmoscoll(self, wrsbox, seasonindex1, seasonindex2, qltyband, quantile):
        self.collkwargs.update(filterpoly=wrsbox.centroid().geometry())
        seascoll = Collection(collkwargs=self.collkwargs, **self.seasargs)
        seascoll.update_basecoll(self.updater)
        clipper = wrsbox.geometry().buffer(-200)
        mos = seascoll.biseasonal_quantmos(seasonindex1, seasonindex2, qltyband, quantile).clip(clipper)

        return mos

    def biseasonal_quantmoscoll(self, seasonindex1, seasonindex2, qltyband, quantile):
        _get_mos = lambda box: self._biseasonal_quantmoscoll(ee.Feature(box), seasonindex1, seasonindex2, qltyband,
                                                             quantile)
        moscoll = self.wrsboxes.map(_get_mos)
        filt = ee.Filter.And(ee.Filter.neq('bands_1', 0), ee.Filter.neq('bands_2', 0))
        # dropping composites without full set of bands
        return moscoll.filter(filt)

"""
Port from users/georgeazzari/EEtools:seasonal.js
"""

s2data = require('users/georgeazzari/EEtools:s2.data.js')
s1tools = require('users/georgeazzari/EEtools:s1.data.tools.js')
fsetrees = require('users/georgeazzari/EEtools:s2.cloudtree.fse.africa.js')


def getS1Plus(poly, year, correctlia, addspeckle):
    startdate = ee.Date.fromYMD(year - 1, 10, 1)
    enddate = ee.Date.fromYMD(year, 10, 1)
    # correctlia = True
    addbands = True
    # addspeckle = True
    addtexture = False
    orbit = None
    s1 = s1tools.getS1IWH(poly, startdate, enddate, correctlia, addbands, addspeckle, addtexture, orbit)
    # Filter to get images with VV and VH dual polarization.
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .map(imgtools.addDOY)
        .map(function(img){
            polbands = s1tools.getPolBands(img)
            img = img.addBands(s1tools.toNatural(img.select(polbands)), null, true)
            return img
        })

    # NOTE: the output of this had polarization bands in NATURAL units
    return s1





"""/*------------------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------- Get S2 collection ---------------------------------------------------*/"""


def getS2Plus(region, year, addvis):
    startdate = ee.Date.fromYMD(year - 1, 10, 1)
    enddate = ee.Date.fromYMD(year, 10, 1)

    s2plus = s2data.getS2(startdate, enddate, region, false, false, false, true).map(

        function(img){

            # Masking with first version of our decision tree
            fsemask = fsetrees.decisionTreeclass(img).select([0],['QA_FSEV1']).eq(1)
            img  = img.updateMask(fsemask).select(s2data.opticalbands)

            # Rescaling to 0-1
            img = s2data.rescaleS2(img)

            return img

        })
    
    if(addvis===true){
      # Adding VIs
      
      s2plus = s2plus.map(
        function(img){
      
          img = s2data.addSWVIs(img)
          img = s2data.addRededgeExtras(img)
          
          return img
          
      })
    }

    return s2plus


"""/*---------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------- TEMPORAL COMPOSITES -----------------------------------------*/"""


def getNewNames(bandnames, suffix):
    newnames = bandnames.map(
        lambda val: ee.String(val)
        .cat(ee.String("_"))
        .cat(ee.String(suffix))
    )
    return newnames


def renameBands(img, suffix):
    bandnames = img.bandNames()
    newnames = getNewNames(bandnames, suffix)
    return img.select(bandnames, newnames)


def appendSeasonBand(current, previous):
    # Rename the band
    current = renameBands(current, current.get('season'))
    # Append it to the result (Note: only return current item on first element/iteration)
    accum = ee.Algorithms.If(ee.Algorithms.IsEqual(previous, None), current, current.addBands(ee.Image(previous)))
    # Return the accumulation
    return accum


#  /* Compute seasonal median composites for a generic collection */
def seasonMedians(imgcoll, bandnames, year, asimgcoll, addtexture, glcmvars):

    f0a = ee.Filter.and(ee.Filter.inList('MONTH', [9, 10, 11]), ee.Filter.eq('YEAR', year - 1))  # Oct, Nov, Dec
    f0b = ee.Filter.and(ee.Filter.inList('MONTH', [0]), ee.Filter.eq('YEAR', year))  # Jan
    f1 = ee.Filter.and(ee.Filter.inList('MONTH', [1, 2, 3, 4]), ee.Filter.eq('YEAR', year))
    f2 = ee.Filter.and(ee.Filter.inList('MONTH', [5, 6, 7, 8]), ee.Filter.eq('YEAR', year))

    bim0 = imgcoll.select(bandnames).filter(ee.Filter.or(f0a, f0b)).median()
    bim1 = imgcoll.select(bandnames).filter(f1).median()
    bim2 = imgcoll.select(bandnames).filter(f2).median()

    if asimgcoll:
        return ee.ImageCollection.fromImages([
            bim0.set({
                'nbands': bim0.bandNames().size(),
                'season': 'S1',
                # The DOY band is quite useless for seasonal composites and complicates things afterwords.
                'BANDNAMES': getNewNames(bim0.bandNames().removeAll(['DOY']), 'S1'),
                'DOYNAMES': getNewNames(['DOY'], 'S1')  # This is a dummy for glcm compositing
            }),
            bim1.set({
                'nbands': bim1.bandNames().size(),
                'season': 'S2',
                # The DOY band is quite useless for seasonal composites and complicates things afterwords.
                'BANDNAMES': getNewNames(bim1.bandNames().removeAll(['DOY']), 'S2'),
                'DOYNAMES': getNewNames(['DOY'], 'S2')  # This is a dummy for glcm compositing
            }),
            bim2.set({
                'nbands': bim2.bandNames().size(),
                'season': 'S3',
                # The DOY band is quite useless for seasonal composites and complicates things afterwords.
                'BANDNAMES': getNewNames(bim2.bandNames().removeAll(['DOY']), 'S3'),
                'DOYNAMES': getNewNames(['DOY'], 'S3')  # This is a dummy for glcm compositing
            })
        ]).filter(ee.Filter.neq('nbands', 0))
    else:
        medcomp = renameBands(bim0, 'S1') \
            .addBands(renameBands(bim1, 'S2')) \
            .addBands(renameBands(bim2, 'S3')) \
            .set('year', year)

        if addtexture:
            scaler = ee.Dictionary.fromLists(
                medcomp.bandNames(),
                ee.List.repeat(1e4, medcomp.bandNames().length())
            )
            glcm = imgtools.getGLCMTexture(medcomp, 4, None, True, scaler)

            if glcmvars is not None:
                selbands = ee.List(glcmvars).map(
                    lambda n: return getNewNames(medcomp.bandNames(), n)
                ).flatten()
                glcm = glcm.select(selbands)

            medcomp = medcomp.addBands(glcm)

        return medcomp


def s1Medians(region, year, correctlia, addspeckle, addtexture, glcmvars):
    s1coll = getS1Plus(region, year, correctlia, addspeckle)
    s1bnames = ee.Image(s1coll.first()).bandNames()
    s1med = seasonMedians(s1coll, s1bnames, year, False, addtexture, glcmvars)

    return s1med.set({
        'year': year,
        's1bands': s1bnames,
    })


def s2Medians(region, year, addvis, addtexture, glcmvars):
    s2coll = getS2Plus(region, year, addvis)
    s2bnames = ee.Image(s2coll.first()).bandNames()
    s2med = seasonMedians(s2coll, s2bnames, year, False, addtexture, glcmvars)

    return s2med.set({
        'year': year,
        's2bands': s2bnames
    })


def sentinel_combined_medians(region, year, addvis, correctlia, addspeckle, addtexture, glcmvars):
    s1med = s1Medians(region, year, correctlia, addspeckle, addtexture, glcmvars)
    s2med = s2Medians(region, year, addvis, addtexture, glcmvars)
    return ee.Image(s2med.addBands(s1med).copyProperties(s1med))
