"""
Author: Anthony Perez

Image datasources related to human settlement patterns (population, urbanization, etc.)
"""
import ee
from gee_tools.datasources.interface import MultiImageDatasource, GlobalImageDatasource, SingleImageDatasource, DatasourceError


def get_closest_to_date(img_coll, date):

    def closest(current_scene, closest_scene):
        current_scene = ee.Image(current_scene)
        closest_scene = ee.Image(closest_scene)
        current_date = ee.Date(current_scene.date())
        closest_date = ee.Date(closest_scene.date())
        curr_diff = ee.Number(date.difference(current_date, 'day')).abs()
        closest_diff = ee.Number(date.difference(closest_date, 'day')).abs()
        return ee.Algorithms.If(closest_diff.gte(curr_diff),
                                current_scene,
                                closest_scene)

    return ee.Image(img_coll.iterate(closest, img_coll.first()))


def get_center_date(start_date, end_date):
    start_date = ee.Date(start_date)
    end_date = ee.Date(end_date)
    diff = ee.Number(end_date.difference(start_date, 'day'))
    advance = ee.Number(diff).divide(ee.Number(2.0))
    return start_date.advance(advance, 'day')


class GHSLPop(GlobalImageDatasource):
    """
    Image Properties

    POPULATION    Number of people per cell

    kwargs:
    :param year:  Sets the start and end date to the start and end of the given year
    :param use_closest_image:  Selects the image closest to the center of the start and end date range
    """

    def build_img_coll(self, year=None, use_closest_image=False):
        if year is not None:
            if not isinstance(year, int):
                raise ValueError('Expected year to be integral.  Got: {}'.format(type(year)))
            self.start_date = str(year) + '-1-1'
            self.end_date = str(year) + '-12-31'

        self.pop = ee.ImageCollection('JRC/GHSL/P2016/POP_GPW_GLOBE_V1')

        if use_closest_image:
            center_date = get_center_date(self.start_date, self.end_date)
            self.pop = ee.ImageCollection(get_closest_to_date(self.pop, center_date))
        else:
            self.pop = self.pop.filterDate(self.start_date, self.end_date)

        self.pop = self.pop.map(self.rename_GHLS_pop) \
                       .sort('system:time_start')

    def get_img_coll(self):
        return self.pop

    @staticmethod
    def rename_GHLS_pop(scene):
        band_names = ['POPULATION']
        return scene.select(list(range(len(band_names))), band_names)


class GHSLUrban(GlobalImageDatasource):
    """
    Image Properties

    SMOD	    Degree of urbanization (See table below).

    Value	Color	Description
    0	    000000	Inhabited areas
    1	    448564	RUR (rural grid cells)
    2	    70daa4	LDC (low density clusters)
    3	    ffffff	HDC (high density clusters)

    If separate is set to True, the resulting images have separate bands
    for each value above.  Each pixel has a value of 1 iff the corresponding pixel in the
    original band matches the new band's corresponding value.

    If separate, band names are ["INHABITED", "RUR", "LDC", "HDC"]
    """

    def build_img_coll(self, year=None, separate=False, use_closest_image=False):
        if year is not None:
            if not isinstance(year, int):
                raise ValueError('Expected year to be integral.  Got: {}'.format(type(year)))
            self.start_date = str(year) + '-1-1'
            self.end_date = str(year) + '-12-31'

        self.urban = ee.ImageCollection('JRC/GHSL/P2016/SMOD_POP_GLOBE_V1')

        if use_closest_image:
            center_date = get_center_date(self.start_date, self.end_date)
            self.urban = ee.ImageCollection(get_closest_to_date(self.urban, center_date))
        else:
            self.urban = self.urban.filterDate(self.start_date, self.end_date)

        if separate:
            self.urban = self.urban.map(self.separate_urban_bands)
        else:
            self.urban = self.urban.map(self.rename_GHLS_Urban)

        self.urban = self.urban.sort('system:time_start')

    def get_img_coll(self):
        return self.urban

    @staticmethod
    def rename_GHLS_Urban(scene):
        band_names = ["SMOD"]
        return scene.select(list(range(len(band_names))), band_names)

    @staticmethod
    def separate_urban_bands(scene):
        """
        Value	Color	Description
        0	    000000	Inhabited areas
        1	    448564	RUR (rural grid cells)
        2	    70daa4	LDC (low density clusters)
        3	    ffffff	HDC (high density clusters)
        """
        smod = scene.select([0])

        inhabited = smod.eq(0).select([0], ["INHABITED"])
        rural = smod.eq(1).select([0], ["RUR"])
        low_denisty = smod.eq(2).select([0], ["LDC"])
        high_denisty = smod.eq(3).select([0], ["HDC"])

        split_bands = scene.addBands([inhabited, rural, low_denisty, high_denisty])
        split_bands = split_bands.select(list(range(1, 5)))
        return split_bands

    @staticmethod
    def get_band_names(separate=False):
        """Return a list of the output band names."""
        if separate:
            return ["INHABITED", "RUR", "LDC", "HDC"]
        return ['SMOD']


class CityAccessibility(SingleImageDatasource):
    """Oxford/MAP/accessibility_to_cities_2015_v1_0"""

    def build_img_coll(self):
        self.dist_to_road = ee.Image("Oxford/MAP/accessibility_to_cities_2015_v1_0")
        self.dist_to_road = self.dist_to_road.select(['accessibility'], ['ACCESSIBILITY'])
        self.coll = ee.ImageCollection([self.dist_to_road])

    def get_img_coll(self):
        return self.coll
