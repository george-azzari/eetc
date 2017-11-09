import ee
import arraycomposites
import lndsatimgtools as imgtools


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
