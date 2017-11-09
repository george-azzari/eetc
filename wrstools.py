import ee
import arraycomposites


def split_wrscode(f):
    wrscode = ee.String(f.get('name'))
    wrs_spl = wrscode.split('_')
    p = ee.String(wrs_spl.get(0))
    r = ee.String(wrs_spl.get(1))
    wrslab = ee.String('P').cat(p).cat(ee.String('R').cat(r))
    # id = p.multiply(10000).add(r)
    return f.set({'WRS_PATH': ee.Number.parse(p),
                  'WRS_ROW': ee.Number.parse(r),
                  'WRSLAB': wrslab})


wrs = ee.FeatureCollection('ft:1VQQdn0RisuaRdQ7ezcenk7rzD9ff6icUYJEgG1GZ').map(split_wrscode)


def fetch_wrsboxes(wrscoll, poly):
    return wrscoll.filterBounds(poly).map(split_wrscode)


def groupby_row(wrscoll):
    distinct_paths = wrscoll.distinct("WRS_ROW")
    joined_rows = ee.Join.saveAll("matches").apply(distinct_paths, wrscoll, ee.Filter.equals("WRS_ROW", None, "WRS_ROW"))
    return joined_rows


def groupby_path(wrscoll):
    distinct_paths = wrscoll.distinct("WRS_PATH")
    joined_rows = ee.Join.saveAll("matches").apply(distinct_paths, wrscoll, ee.Filter.equals("WRS_PATH", None, "WRS_PATH"))
    return joined_rows  # rows are stored in the property 'matches' of each feature (as ee.List)


def filterby_pathrow(p, r, wrsfetcoll):
    p = ee.Number(p)
    r = ee.Number(r)
    f = ee.Filter.And(ee.Filter.eq('WRS_ROW', r), ee.Filter.eq('WRS_PATH', p))
    return wrsfetcoll.filter(f)


def get_left(wrsfeat, wrsfetcoll):
    # Note: explicit casting of wrsfeat to a specific object becomes necessary when this method gets mapped (e.g. from
    #       the iterate_across methods). Notice that if this method was simply iterated (e.g. from the iterate_along
    #       methods) the explicit cast would not be necessary. However, it does not seem to make a difference whether
    #       wrsfeat is casted to ee.Feature or ee.Image (e.g. the mean-correction method works both ways) so I will
    #       cast it to a feature, as it seems more generic to me.
    wrsfeat = ee.Feature(wrsfeat)
    r = ee.Number(wrsfeat.get('WRS_ROW'))
    p = ee.Number(wrsfeat.get('WRS_PATH'))
    sel = wrsfetcoll.filter(ee.Filter.And(ee.Filter.eq('WRS_ROW', r),
                                          ee.Filter.gte('WRS_PATH', p))).sort('WRS_PATH', True)
    return sel.set({'WRS_ROW': r})


def get_right(wrsfeat, wrsfetcoll):
    # Note: explicit casting of wrsfeat to a specific object becomes necessary when this method gets mapped (e.g. from
    #       the iterate_across methods). Notice that if this method was simply iterated (e.g. from the iterate_along
    #       methods) the explicit cast would not be necessary. However, it does not seem to make a difference whether
    #       wrsfeat is casted to ee.Feature or ee.Image (e.g. the mean-correction method works both ways) so I will
    #       cast it to a feature, as it seems more generic to me.
    wrsfeat = ee.Feature(wrsfeat)
    r = ee.Number(wrsfeat.get('WRS_ROW'))
    p = ee.Number(wrsfeat.get('WRS_PATH'))
    sel = wrsfetcoll.filter(ee.Filter.And(ee.Filter.eq('WRS_ROW', r),
                                          ee.Filter.lte('WRS_PATH', p))).sort('WRS_PATH', False)
    return sel.set({'WRS_ROW': r})


def get_above(wrsfeat, wrsfetcoll):
    # Note: explicit casting of wrsfeat to a specific object becomes necessary when this method gets mapped (e.g. from
    #       the iterate_across methods). Notice that if this method was simply iterated (e.g. from the iterate_along
    #       methods) the explicit cast would not be necessary. However, it does not seem to make a difference whether
    #       wrsfeat is casted to ee.Feature or ee.Image (e.g. the mean-correction method works both ways) so I will
    #       cast it to a feature, as it seems more generic to me.
    wrsfeat = ee.Feature(wrsfeat)
    r = ee.Number(wrsfeat.get('WRS_ROW'))
    p = ee.Number(wrsfeat.get('WRS_PATH'))
    sel = wrsfetcoll.filter(ee.Filter.And(ee.Filter.eq('WRS_PATH', p),
                                          ee.Filter.lte('WRS_ROW', r))).sort('WRS_ROW', False)
    return sel.set({'WRS_PATH': p})


def get_below(wrsfeat, wrsfetcoll):
    # Note: explicit casting of wrsfeat to a specific object becomes necessary when this method gets mapped (e.g. from
    #       the iterate_across methods). Notice that if this method was simply iterated (e.g. from the iterate_along
    #       methods) the explicit cast would not be necessary. However, it does not seem to make a difference whether
    #       wrsfeat is casted to ee.Feature or ee.Image (e.g. the mean-correction method works both ways) so I will
    #       cast it to a feature, as it seems more generic to me.
    wrsfeat = ee.Feature(wrsfeat)
    r = ee.Number(wrsfeat.get('WRS_ROW'))
    p = ee.Number(wrsfeat.get('WRS_PATH'))
    sel = wrsfetcoll.filter(ee.Filter.And(ee.Filter.eq('WRS_PATH', p),
                                          ee.Filter.gte('WRS_ROW', r))).sort('WRS_ROW', True)
    return sel.set({'WRS_PATH': p})


# /* Given a center box, iterate mean correction along boxes on the left 
#    and the right. Maximum boxes set to 20. Note: can be adapted fot other
#    types of corrections. */
def iterate_alongrow(cntr, wrsboxes, toiterate):
    """
    Given a center box, iterate function toiterate along boxes on the left
    and the right. Maximum boxes set to 50. Note: can be adapted fot other
    types of corrections.
    :param cntr:
    :param wrsboxes:
    :param toiterate:
    :return: ee.List
    """
    start = ee.List([cntr])
    leftlist = ee.List(get_left(cntr, wrsboxes).toList(150)).slice(1)  # slice(1) to drop center box
    rightlist = ee.List(get_right(cntr, wrsboxes).toList(150)).slice(1)  # slice(1) to drop center box
    leftcorrfc = ee.List(leftlist.iterate(toiterate, start))
    rightcorrfc = ee.List(rightlist.iterate(toiterate, start))
    return leftcorrfc.cat(rightcorrfc)  # this is a list


# /* Given a center box, iterate mean correction along boxes above and below.
#    Maximum boxes set to 50. Note: can be adapted fot other types of corrections. */
def iterate_alongpath(cntr, wrsboxes, toiterate):
    """

    :param cntr: a starting feature or image
    :param wrsboxes: a ee.Collection with WRS metadata (either ImageCollection or FeatureCollection)
    :param toiterate: iterable function
    :return: ee.List of features or images
    """
    start = ee.List([cntr])
    abovelist = ee.List(get_above(cntr, wrsboxes).toList(150)).slice(1)  # slice(1) to drop center box
    belowlist = ee.List(get_below(cntr, wrsboxes).toList(150)).slice(1)  # slice(1) to drop center box
    abovecorrfc = ee.List(abovelist.iterate(toiterate, start))
    belowcorrfc = ee.List(belowlist.iterate(toiterate, start))
    return abovecorrfc.cat(belowcorrfc)  # this is a list


def find_pathcenter(wrspathcoll):
    """

    :param wrspathcoll: ee.FeatureCollection/ee.ImageCollection
    :return: ee.Feature/ee.Image
    """
    wrspathcoll = ee.FeatureCollection(wrspathcoll)  # casting won't be necessary if I simply nest in _filter_and_set
    med = ee.Number(wrspathcoll.reduceColumns(ee.Reducer.median(), ee.List(['WRS_ROW'])).get('median')).toInt()
    medrow = ee.Feature(wrspathcoll.filter(ee.Filter.eq('WRS_ROW', med)).first())
    return medrow


# def find_pathcentroid(wrspathcoll):
#     wrspathcoll = ee.FeatureCollection(wrspathcoll)  # casting won't be necessary if I simply nest in _filter_and_set
#     pcntr = ee.Feature(wrspathcoll.union().first()).centroid()
#     wrcntr = ee.Feature(wrspathcoll.filterBounds(pcntr.geometry()).first())
#     return wrcntr


def iterate_alongpathgroup(wrspathgroup, toiterate):
    wrspathcoll = ee.FeatureCollection(ee.List(ee.Feature(wrspathgroup).get('matches')))
    center = ee.Feature(find_pathcenter(wrspathcoll))
    path_iteratedcoll = ee.ImageCollection.fromImages(iterate_alongpath(center, wrspathcoll, toiterate))
    path_mean = path_iteratedcoll.mean().clip(wrspathcoll.union().geometry())
    return path_mean.set({'WRS_PATH': center.get('WRS_PATH'), 'WRS_ROW': 1})


def iterate_alongpathgroups(wrscoll, toiterate):
    pathgroups = groupby_path(wrscoll)
    f = lambda group: iterate_alongpathgroup(group, toiterate)
    return ee.ImageCollection(pathgroups.map(f))


def iterate_acrosspathgroups(wrscoll, pcntr, toiterate):
    pathcoll = iterate_alongpathgroups(wrscoll, toiterate)
    start = ee.Feature(filterby_pathrow(pcntr, 1, pathcoll).first())
    return iterate_alongrow(start, pathcoll, toiterate)


def iterate_acrosspaths(alongpathlist, wrsboxes, toiterate):
    f = lambda b: iterate_alongrow(b, wrsboxes, toiterate)
    lateralpaths = alongpathlist.map(f)
    return lateralpaths.flatten()


def iterate_acrossrows(alongrowlist, wrsboxes, toiterate):
    f = lambda b: iterate_alongpath(b, wrsboxes, toiterate)
    lateralpaths = alongrowlist.map(f)
    return lateralpaths.flatten()


def iterate_acrosswrs(wrsboxes, pcntr, rcntr, alongrows, toiterate):
    wrsboxes_cent = ee.Feature(filterby_pathrow(pcntr, rcntr, wrsboxes).first())
    if alongrows:
        mainbranch = iterate_alongpath(wrsboxes_cent, wrsboxes, toiterate)
        whole = iterate_acrosspaths(mainbranch, wrsboxes, toiterate)
    else:
        mainbranch = iterate_alongrow(wrsboxes_cent, wrsboxes, toiterate)
        whole = iterate_acrossrows(mainbranch, wrsboxes, toiterate)
    return whole


def _get_mean(img, interspoly):
    mean_args = {'reducer': ee.Reducer.mean(),
                 'scale': 240,
                 'tileScale': 8,
                 'bestEffort': False,
                 'geometry': interspoly,
                 'maxPixels': 1e13}
    return ee.Number(img.select([0], ['mean']).reduceRegion(**mean_args).get('mean'))


def _get_mean_corrfactor(refimg, destimg, interspoly):
    refmean = _get_mean(ee.Image(refimg), interspoly)
    destmean = _get_mean(ee.Image(destimg), interspoly)
    deltamean = ee.Algorithms.If(refmean, ee.Algorithms.If(destmean, refmean.subtract(destmean), 0), 0)
    return ee.Number(deltamean)


def _correct_bymean(refimg, destimg, interspoly):
    deltamean = _get_mean_corrfactor(refimg, destimg, interspoly)
    corrimg = destimg.add(deltamean).copyProperties(destimg)
    # NOTE: resulting image looses footprint info (no system:footprint" metadata); had to add a clip.
    #       This would not be necessary with the WRS-box approach.
    return ee.Image(corrimg).clip(destimg.geometry())


def _iterable_meancorr(wrsimg, imglist):
    imglist = ee.List(imglist)
    wrsimg = ee.Image(wrsimg)
    refimg = ee.Image(imglist.get(-1))
    # Note: here is one point in which the question I raised in method "biseasonal_wrstilescym" will come
    #       into play. If the scym image was clipped, I can grab directly its geometry, otherwise I will
    #       have to call the corresponding WRS box.
    # ------------------------------------------------------------------------
    # Here is the flow assuming that intersection is found using corresponding WRS boxes:
    # destpoly = wrstools.filterby_pathrow(ee.Number(wrsimg.get('WRS_PATH')), ee.Number(wrsimg.get('WRS_ROW')),
    #                                      wrstools.wrs)
    # refpoly = wrstools.filterby_pathrow(ee.Number(refimg.get('WRS_PATH')), ee.Number(refimg.get('WRS_ROW')),
    #                                     wrstools.wrs)
    # dest_r = ee.Number(wrsimg.get('WRS_ROW'))
    # dest_p = ee.Number(wrsimg.get('WRS_PATH'))
    # ------------------------------------------------------------------------
    # Here is the flow assuming wrsimg is pre-clipped and has a WRS geometry:
    destpoly = wrsimg.geometry()
    refpoly = refimg.geometry()
    corrimg = _correct_bymean(refimg, wrsimg, refpoly.intersection(destpoly))
    return imglist.add(ee.Image(corrimg))


def _correct_wrsgrid_bymean(wrscoll, pcntr, rcntr, alongrows):
    corrcoll = iterate_acrosswrs(wrscoll, pcntr, rcntr, alongrows, _iterable_meancorr)
    return corrcoll


def _blend_diff(wrsimg, diff):
    blendedimg = ee.Image(wrsimg.select('yield').subtract(diff)).toFloat()
    return wrsimg.addBands(blendedimg, ee.List(['yield']), True)


def _blend(wrsimg, modimg):
    wrsimg = ee.Image(wrsimg)
    modimg = ee.Image(modimg)
    diffimg = ee.Image(wrsimg.select('yield').subtract(modimg.select('yield'))).select([0], ['Dyield'])
    reduced = diffimg.reduceRegion(reducer=ee.Reducer.median(), scale=960, bestEffort=False)
    diff = ee.Number(ee.Dictionary(reduced).get('Dyield'))
    # TODO: I did'n find a way around this If yet. Perhaps with pre-mapping diffimg across a collection with dropnulls?
    corr = ee.Algorithms.If(diff, _blend_diff(wrsimg, diff), None)
    # return wrsimg.addBands(ee.Image(corr), ee.List(['yield']), True)
    # TODO: should probably add a toInt16() to the output as offsetting may convert it to double.
    return ee.Image(corr)


def modis_blending_byyear(modiscoll, wrscoll, year):
    modimg = ee.Image(modiscoll.filter(ee.Filter.eq('year', year)).first())
    blend = lambda wrsimg: _blend(wrsimg, modimg)
    # TODO: I could not find a way around this toList. Should ask the devs.
    corrcoll = ee.ImageCollection.fromImages(wrscoll.map(blend, True).toList(150))
    return arraycomposites.nthquantile_mosaic(corrcoll, 'yield', 0.5).set({'year': year})
    # return ee.Image(corrcoll.median()).set({'year': year})
