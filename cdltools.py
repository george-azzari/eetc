"""
Author: George Azzari (gazzari@stanford.edu)
Center on Food Security and the Environment
Department of Earth System Science
Stanford University
"""


import ee
import imgtools
import export_tables as tabs

cdlcoll = ee.ImageCollection("USDA/NASS/CDL")

states = ee.FeatureCollection("ft:17aT9Ud-YnGiXdXEJUyycH2ocUqreOeKGbzCkUw")
counties = ee.FeatureCollection("ft:18Ayj5e7JxxtTPm1BdMnnzWbZMrxMB49eqGDTsaSp")


def get_cdl(year):

    cdlcoll = ee.ImageCollection("USDA/NASS/CDL")

    cdl = ee.Image(cdlcoll.filter(ee.Filter.calendarRange(year, year, 'YEAR')).first())
    cdl = cdl.select(['cropland'], ['CDL'])

    cropmask = cdl.gte(196).add(cdl.lte(60))
    cropmask = cropmask.add(cdl.gte(66).multiply(cdl.lte(77))).select([0], ['CROPMASK'])

    tosample = imgtools.add_latlon(cdl.updateMask(cropmask).addBands(cropmask))

    return tosample


def format_rand_id(sample, suffix):

    r = ee.String(ee.Number(sample.get('rand')).multiply(1000000).toInt().format())
    f = ee.String(ee.Number(sample.get('FIPS formula')).toInt().format())
    s = ee.String(sample.get('State Abbr'))

    pointid = s.cat(ee.String('-')).cat(f).cat(ee.String('-')).cat(r).cat(ee.String(suffix))

    return sample.set('SAMPLEID', pointid)


def _get_cdl_info(cdlpoint, vals, names, colors):

    tval = ee.Number(cdlpoint.get('CDL'))
    valindx = vals.indexOf(tval)

    return cdlpoint.set('CDL_NAME', ee.String(names.get(valindx)),
                        'CDL_COLOR', ee.String(colors.get(valindx)))


def get_cdl_stratgrid(region, year, npoints, random_seed):

    cdl = get_cdl(year)

    samples = cdl.stratifiedSample(
        numPoints=npoints,
        classBand='CROPMASK',
        region=region.geometry(),
        scale=30,
        seed=random_seed,
        classValues=[0, 1],
        classPoints=[0, npoints],
        tileScale=16)

    # Copy properties from region feature and add year
    samples = samples.map(lambda s: ee.Feature(s.copyProperties(region)).set('year', year))

    # Add class name and color based on value:
    vals = ee.List(cdl.get('cropland_class_values'))
    names = ee.List(cdl.get('cropland_class_names'))
    colors = ee.List(cdl.get('cropland_class_palette'))
    samples = samples.map(lambda p: _get_cdl_info(p, vals, names, colors))

    # Finally, add unique ID based on county fips, random number, and year
    samples = samples.randomColumn('rand', 3718).map(lambda p: format_rand_id(p, ee.String('-'+str(year))))

    return samples


def export_cdl_stratgrid(features, year, npoints, random_seed, fname, export_to='drive'):

    grids = features.map(lambda county: get_cdl_stratgrid(county, year, npoints, random_seed)).flatten()

    task = tabs.export_features(grids, fname, export_to)

    return task


"""
EXAMPLE:

selstates = ['ND', 'SD', 'NE', 'KS', 'MO', 'IA', 'MN', 'WI', 'IL', 'IN', 'MI','OH', 'KY']
selcounties = counties.filter(ee.Filter.inList('State Abbr', selstates))

fname = 'cdlgrid_midwest_stratbycounty_500pts_r345_2016'
grid_bycounty = cdltools.export_cdl_stratgrid(selcounties, 2016, 500, 345, fname, export_to='drive')

"""